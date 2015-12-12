class Camera;
class Camera
{
    public:
        float azimuth, radius, zenith;
        Vec3D eye, eyeVector, pos;
        Mat4 view;
};

class Mat4;
class Mat4
{
    float mat[4][4];
public:
    Mat4(Point3D p1, Point3D p2, Point3D p3, Point3D p4);
    Mat4(Vec3D e, Vec3D v, Vec3D u);

};
class Vec3D;
class Vec3D
{
public:
    float x, y, z;
    Vec3D();
    Vec3D(float xx, float yy, float zz);
    float dot(Vec3D rhs); // dot produkt
    Vec3D mul(float d); // nasobeni skalarem
    Vec3D normalized(); // normalizace
    Vec3D cross(Vec3D v); // cross produkt
};

class Point3D;
class Point3D
{
public:
    float x, y, z, w;
    Point3D(float xx, float yy, float zz);

};