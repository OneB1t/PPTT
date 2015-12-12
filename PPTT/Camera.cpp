#include "Camera.h";
#include "PPTT_core.h"
#include <cmath>

Camera::Camera()
{
    azimuth = 5.0f;
    zenith = 0.0f;
    radius = 1.0f;
    eyeVector = Vec3D(0,1,0);
    pos = Vec3D(50.0f, 150.0f, 50.0f);
    computeCameraMatrix();
}

void Camera::computeCameraMatrix()
{
    eyeVector = Vec3D(cos(azimuth) * cos(zenith),sin(azimuth) * cos(zenith),sin(zenith));
    from = Vec3D(pos);
    to = eyeVector.mul(radius);
    upV.x = (cos(azimuth) * cos(zenith + PI / 2));
    upV.y = sin(azimuth) * cos(zenith + PI / 2);
    upV.z = sin(zenith + PI / 2);
}

void Camera::addAzimuth(float ang)
{
    azimuth += ang;
    computeCameraMatrix();
}

void Camera::addRadius(float dist)
{
    if(radius + dist < 0.1f)
        return;
    radius += dist;
    computeCameraMatrix();
}

void Camera::addZenith(float ang)
{
    if(abs(zenith + ang) <= PI / 2) {
        zenith += ang;
        computeCameraMatrix();
    }
}

void Camera::backward(float speed)
{
    forward((-1) * speed);
}

void Camera::forward(float speed)
{
    pos = pos.add(Vec3D(cos(azimuth) *cos(zenith), sin(azimuth) * cos(zenith), sin(zenith) * speed)); // maybe issue here need to check
    computeCameraMatrix();
}

void Camera::left(float speed)
{
    right((-1) * speed);
}

void Camera::right(float speed)
{
    pos = pos.add(Vec3D(cos(azimuth - PI / 2),sin(azimuth - PI / 2), 0.0f).mul(speed));
    computeCameraMatrix();
}

void Camera::down(float speed)
{
    pos.z -= speed;
    computeCameraMatrix();
}

void Camera::up(float speed)
{
    pos.z += speed;
    computeCameraMatrix();
}
void Camera::move(Vec3D dir)
{
    pos = pos.add(dir);
    computeCameraMatrix();
}
void Camera::moveback(Vec3D dir)
{
    dir.x = -dir.x;
    dir.y = -dir.y;
    dir.z = -dir.z;
    pos = pos.add(dir);
    computeCameraMatrix();
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

Vec3D Vec3D::add(Vec3D v)
{
    return Vec3D(x + v.x, y + v.y, z + v.z);
}

Point3D::Point3D(float xx, float yy, float zz)
{
    x = xx;
    y = yy;
    z = zz;
    w = 1.0f;
}
