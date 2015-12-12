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
    Vec3D add(Vec3D v);
};

class Camera;
class Camera
{
    public:
        float azimuth, radius, zenith;
        int ox,oy;
        Vec3D eye;
        Vec3D from,to,upV;
        Vec3D eyeVector;
        Vec3D pos;
        Camera();
        void computeCameraMatrix();
        void addAzimuth(float ang);
        void addRadius(float dist);
        void addZenith(float ang);
        void backward(float speed);
        void forward(float speed);
        void left(float speed);
        void right(float speed);
        void down(float speed);
        void up(float speed);
        void move(Vec3D dir);
        void moveback(Vec3D dir);
};

class Point3D;
class Point3D
{
public:
    float x, y, z, w;
    Point3D(float xx, float yy, float zz);

};