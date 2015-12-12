#include "Camera.h";
#include <math.h>;

Mat4::Mat4(Point3D p1, Point3D p2, Point3D p3, Point3D p4)
{
    mat[0][0] = p1.x;
    mat[0][1] = p1.y;
    mat[0][2] = p1.z;
    mat[0][3] = p1.w;
    mat[1][0] = p2.x;
    mat[1][1] = p2.y;
    mat[1][2] = p2.z;
    mat[1][3] = p2.w;
    mat[2][0] = p3.x;
    mat[2][1] = p3.y;
    mat[2][2] = p3.z;
    mat[2][3] = p3.w;
    mat[3][0] = p4.x;
    mat[3][1] = p4.y;
    mat[3][2] = p4.z;
    mat[3][3] = p4.w;
}

Mat4::Mat4(Vec3D e, Vec3D v, Vec3D u)
{
    Vec3D x, y, z;
    z = v.mul(-1.0f).normalized();
    x = u.cross(z).normalized();
    y = z.cross(x);
    mat[0][0] = x.x;
    mat[1][0] = x.y;
    mat[2][0] = x.z;
    mat[3][0] = -e.dot(x);
    mat[0][1] = y.x;
    mat[1][1] = y.y;
    mat[2][1] = y.z;
    mat[3][1] = -e.dot(y);
    mat[0][2] = z.x;
    mat[1][2] = z.y;
    mat[2][2] = z.z;
    mat[3][2] = -e.dot(z);
}

Vec3D::Vec3D(float xx, float yy, float zz)
{
    x = xx;
    y = yy;
    z = zz;
}

Vec3D::Vec3D()
{
    x = y = z = 0.0f;
}

float Vec3D::dot(Vec3D rhs)
{
    return x * rhs.x + y * rhs.y + z * rhs.z;
}

Vec3D Vec3D::mul(float d)
{
    return Vec3D(x * d, y * d, z * d);
}

Vec3D Vec3D::normalized()
{
    double len = sqrt(x*x + y*y + z*z);
    if(len == 0.0f)
        return Vec3D(0, 0, 0);
    return Vec3D(x / len, y / len, z / len);
}

Vec3D Vec3D::cross(Vec3D v)
{
    return Vec3D(y * v.z - z * v.y, z * v.x - x * v.z, x
        * v.y - y * v.x);
}

Point3D::Point3D(float xx, float yy, float zz)
{
    x = xx;
    y = yy;
    z = zz;
    w = 1.0f;
}
