/* I had an idea, once */

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>

typedef uint32_t u32;
typedef float f32;

struct Vec3 {
    f32 x, y, z;
};

Vec3
operator+(const Vec3 &rhs, const Vec3 &lhs) {
    Vec3 v;
    v.x = rhs.x + lhs.x;
    v.y = rhs.y + lhs.y;
    v.z = rhs.z + lhs.z;
    return v;
}

Vec3
operator-(const Vec3 &rhs, const Vec3 &lhs) {
    Vec3 v;
    v.x = rhs.x - lhs.x;
    v.y = rhs.y - lhs.y;
    v.z = rhs.z - lhs.z;
    return v;
}

Vec3
operator*(const Vec3 &rhs, f32 s) {
    Vec3 v;
    v.x = rhs.x * s;
    v.y = rhs.y * s;
    v.z = rhs.z * s;
    return v;
}

f32
dot(const Vec3 &rhs, const Vec3 &lhs) {
    return rhs.x * lhs.x + rhs.y * lhs.y + rhs.z * lhs.z;
}

Vec3
normalize(const Vec3 &v) {
    f32 length = sqrtf(dot(v, v));
    Vec3 n;
    n.x = v.x / length;
    n.y = v.y / length;
    n.z = v.z / length;
    return n;
}

Vec3
cross(const Vec3 &rhs, const Vec3 &lhs) {
    Vec3 c;
    c.x = rhs.y * lhs.z - rhs.z * lhs.y;
    c.y = rhs.z * lhs.x - rhs.x * lhs.z;
    c.z = rhs.x * lhs.y - rhs.y * lhs.x;
    return c; 
}

struct Entity {
    Vec3 position;
    f32 radius;
};

struct WorldRegion {
    Entity *entities;
    size_t  num_entities;
    Vec3 position;
    f32  size; // assume that width,height,depth are equal
    // and that the region is axis aligned
    Vec3 vertices[8];
};

struct Plane {
    // plane equation
    // nx * x + ny * y + nz * z + d = 0
    f32 nx, ny, nz, d;
};

struct Camera {
    Plane frustum[6];
};

size_t
cull_world_regions(const Camera *camera, const WorldRegion *regions,
                    size_t num, WorldRegion *out) {
    u32 vertex_luts[6][8];
    // Precompute the vertex luts for each frustum plane 
    for (int planeidx = 0; planeidx < 6; planeidx++) {
        
    }
    // This assumes that the camera planes are in world space
    size_t outidx = 0;
    Plane *const frustum = (Plane *const)camera->frustum;
    for (size_t i = 0; i < num; i++) {
        bool outside = false;
        // NOTE(Kevin): A decent compiler should unroll this loop
        for (int planeidx = 0; planeidx < 6; planeidx++) {
            // NOT(Kevin): I think this can be done via bit operations.
            u32 signs = (frustum[planeidx].nx < 0.0f) |
                        ((frustum[planeidx].ny < 0.0f) << 1) |
                        ((frustum[planeidx].nz < 0.0f) << 2);
            printf("signs: %f %f %f -> %u\n",
                    frustum[planeidx].nx,
                    frustum[planeidx].ny,
                    frustum[planeidx].nz,
                    signs);
            // negative point
            Vec3 n = regions[i].vertices[vertex_luts[planeidx][signs]];
            f32 a = n.x * frustum[planeidx].nx + n.y * frustum[planeidx].ny + n.z * frustum[planeidx].nz + frustum[planeidx].d;
            outside = outside || (a < 0.0f);
        } 
        out[outidx] = regions[i];
        outidx += outside;
    } 
    return outidx;
}

size_t
cull_entities(const Camera *camera, const WorldRegion *regions, size_t num, Entity *out) {
}

void
update_camera_frustum(Vec3 position, Vec3 up, Vec3 forward, Vec3 right, f32 near, f32 far, f32 fov, f32 aspect, Camera *out) { 
    Vec3 camX = { position.x, 0.0f, 0.0f };
    Vec3 camY = { 0.0f, position.y, 0.0f };
    Vec3 nearCenter = position - forward * near;
    Vec3 farCenter = position - forward * far;
    f32 nearHeight = 2 * tanf(fov/ 2) * near;
    f32 farHeight = 2 * tanf(fov / 2) * far;
    f32 nearWidth = nearHeight * aspect;
    f32 farWidth = farHeight * aspect;
    Vec3 farTopLeft = farCenter + up * (farHeight*0.5f) - right * (farWidth*0.5f);
    Vec3 farTopRight = farCenter + up * (farHeight*0.5f) + right * (farWidth*0.5f);
    Vec3 farBottomLeft = farCenter - up * (farHeight*0.5f) - right * (farWidth*0.5f);
    Vec3 farBottomRight = farCenter - up * (farHeight*0.5f) + right * (farWidth*0.5f);

    Vec3 nearTopLeft = nearCenter + camY * (nearHeight*0.5f) - camX * (nearWidth*0.5f);
    Vec3 nearTopRight = nearCenter + camY * (nearHeight*0.5f) + camX * (nearWidth*0.5f);
    Vec3 nearBottomLeft = nearCenter - camY * (nearHeight*0.5f) - camX * (nearWidth*0.5f);
    Vec3 nearBottomRight = nearCenter - camY * (nearHeight*0.5f) + camX * (nearWidth*0.5f);
    Vec3 p0, p1, p2;

    p0 = nearBottomLeft; p1 = farBottomLeft; p2 = farTopLeft;
    Vec3 leftPlaneNormal = normalize(cross(p1-p0, p2-p1));
    f32 leftPlaneOffset = dot(leftPlaneNormal, p0);

    p0 = nearTopLeft; p1 = farTopLeft; p2 = farTopRight;
    Vec3 topPlaneNormal = normalize(cross(p1-p0, p2-p1));
    f32 topPlaneOffset = dot(topPlaneNormal , p0);

    p0 = nearTopRight; p1 = farTopRight; p2 = farBottomRight;
    Vec3 rightPlaneNormal = normalize(cross(p1-p0, p2-p1));
    f32 rightPlaneOffset = dot(rightPlaneNormal , p0);

    p0 = nearBottomRight; p1 = farBottomRight; p2 = farBottomLeft;
    Vec3 bottomPlaneNormal = normalize(cross(p1-p0, p2-p1));
    f32 bottomPlaneOffset = dot(bottomPlaneNormal , p0);

    Vec3 nearPlaneNormal = forward * -1.0f;
    f32 nearPlaneOffset = dot(nearPlaneNormal, nearCenter);

    Vec3 farPlaneNormal = forward;
    f32 farPlaneOffset  = dot(farPlaneNormal, farCenter);
    
    out->frustum[0].nx = nearPlaneNormal.x;
    out->frustum[0].ny = nearPlaneNormal.y;
    out->frustum[0].nz = nearPlaneNormal.z;
    out->frustum[0].d  = nearPlaneOffset;
    out->frustum[1].nx = farPlaneNormal.x;
    out->frustum[1].ny = farPlaneNormal.y;
    out->frustum[1].nz = farPlaneNormal.z;
    out->frustum[1].d  = farPlaneOffset;
    out->frustum[2].nx = leftPlaneNormal.x;
    out->frustum[2].ny = leftPlaneNormal.y;
    out->frustum[2].nz = leftPlaneNormal.z;
    out->frustum[2].d  = leftPlaneOffset;
    out->frustum[3].nx = topPlaneNormal.x;
    out->frustum[3].ny = topPlaneNormal.y;
    out->frustum[3].nz = topPlaneNormal.z;
    out->frustum[3].d  = topPlaneOffset;
    out->frustum[4].nx = rightPlaneNormal.x;
    out->frustum[4].ny = rightPlaneNormal.y;
    out->frustum[4].nz = rightPlaneNormal.z;
    out->frustum[4].d  = rightPlaneOffset;
    out->frustum[5].nx = bottomPlaneNormal.x;
    out->frustum[5].ny = bottomPlaneNormal.y;
    out->frustum[5].nz = bottomPlaneNormal.z;
    out->frustum[5].d  = bottomPlaneOffset;
}

int
main() {
    WorldRegion regions[1];
    WorldRegion visible_regions[1];
    Camera camera;
    // camera looks in positive z direction
    Vec3 camera_pos = { 0.0f, 0.0f, 0.0f };
    Vec3 camera_up = { 0.0f, 1.0f, 0.0f };
    Vec3 camera_forward = { 0.0f, 0.0f, 1.0f };
    Vec3 camera_right = { 1.0f, 0.0f, 0.0f };
    f32  fov = 45.0f;
    update_camera_frustum(camera_pos, camera_up, camera_forward, camera_right, 0.05f, 1000.0f, fov, 4.0f / 3.0f, &camera);

    Vec3 v = { 0.0f, 0.0f, 0.0f };
    regions[0].position = v; 
    size_t num_visible_regions = cull_world_regions(&camera, regions, 1, visible_regions);

    return 0;
}
