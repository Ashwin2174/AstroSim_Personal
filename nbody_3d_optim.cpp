#include <bits/stdc++.h>
#include <vector>
using namespace std;
const long double G = 6.6743e-11;

class Body {
    public:
        long double mass;
        vector<long double> position;
        vector<long double> momentum;
        long double radius;
        Body(long double m,long double r,long double x,long double y,long double z,long double vx,long double vy,long double vz){
            mass = m;
            radius = r;
            position = {x,y,z};
            momentum = {m*vx,m*vy,m*vz};
        }
};

class solver {
    private:
        void collision(vector<Body>& o, long double e){
            int n = o.size();
            for (int i=0;i<n-1;i++){
                for (int j=i+1;j<n;j++){
                    long double xr = o[j].position[0] - o[i].position[0];
                    long double yr = o[j].position[1] - o[i].position[1];
                    long double d = sqrtl(xr*xr+yr*yr);
                    long double k = o[i].radius+o[j].radius-d;
                    if (d <= (o[i].radius+o[j].radius)){
                        long double delpx = (1+e)*o[i].mass*o[j].mass*(xr*(o[j].momentum[0]/o[j].mass - o[i].momentum[0]/o[i].mass)+yr*(o[j].momentum[1]/o[j].mass - o[i].momentum[1]/o[i].mass))*(xr)/((o[i].mass+o[j].mass)*(d*d));
                        long double delpy = (1+e)*o[i].mass*o[j].mass*(xr*(o[j].momentum[0]/o[j].mass - o[i].momentum[0]/o[i].mass)+yr*(o[j].momentum[1]/o[j].mass - o[i].momentum[1]/o[i].mass))*(yr)/((o[i].mass+o[j].mass)*(d*d));
                        o[i].momentum[0] += delpx;
                        o[i].momentum[1] += delpy;
                        o[j].momentum[0] -= delpx;
                        o[j].momentum[1] -= delpy;
                        o[j].position[0] += k*xr/d;
                        o[j].position[1] += k*yr/d;
                        o[i].position[0] -= k*xr/d;
                        o[i].position[1] -= k*yr/d;
                        cout<<"Collision"<<endl;
                    }
                }
            }
        }   
        vector<long double> px_deriv(vector<Body> o){
            int n = o.size();
            vector<long double> pxd(n,0);
            for (int i=0;i<n-1;i++){
                for (int j=i+1;j<n;j++){
                    long double xr = o[j].position[0] - o[i].position[0];
                    long double yr = o[j].position[1] - o[i].position[1];
                    long double zr = o[j].position[2] - o[i].position[2];
                    pxd[i]+=(G*o[i].mass*o[j].mass*xr)/pow(xr*xr+yr*yr+zr*zr,1.5);
                    pxd[j]-=(G*o[i].mass*o[j].mass*xr)/pow(xr*xr+yr*yr+zr*zr,1.5);
                }
            }
            return pxd;
        }
        vector<long double> py_deriv(vector<Body> o){
            int n = o.size();
            vector<long double> pyd(n,0);
            for (int i=0;i<n-1;i++){
                for (int j=i+1;j<n;j++){
                    long double xr = o[j].position[0] - o[i].position[0];
                    long double yr = o[j].position[1] - o[i].position[1];
                    long double zr = o[j].position[2] - o[i].position[2];
                    pyd[i]+=(G*o[i].mass*o[j].mass*yr)/pow(xr*xr+yr*yr+zr*zr,1.5);
                    pyd[j]-=(G*o[i].mass*o[j].mass*yr)/pow(xr*xr+yr*yr+zr*zr,1.5);
                }
            }
            return pyd;
        }
        vector<long double> pz_deriv(vector<Body> o){
            int n = o.size();
            vector<long double> pzd(n,0);
            for (int i=0;i<n-1;i++){
                for (int j=i+1;j<n;j++){
                    long double xr = o[j].position[0] - o[i].position[0];
                    long double yr = o[j].position[1] - o[i].position[1];
                    long double zr = o[j].position[2] - o[i].position[2];
                    pzd[i]+=(G*o[i].mass*o[j].mass*zr)/pow(xr*xr+yr*yr+zr*zr,1.5);
                    pzd[j]-=(G*o[i].mass*o[j].mass*zr)/pow(xr*xr+yr*yr+zr*zr,1.5);
                }
            }
            return pzd;
        }
        vector<long double> x_deriv(vector<Body> o){
            int n = o.size();
            vector<long double> xd(n,0);
            for (int i=0;i<n;i++){
                xd[i]=o[i].momentum[0]/o[i].mass;
            }
            return xd;
        }
        vector<long double> y_deriv(vector<Body> o){
            int n = o.size();
            vector<long double> yd(n,0);
            for (int i=0;i<n;i++){
                yd[i]=o[i].momentum[1]/o[i].mass;
            }
            return yd;
        }
        vector<long double> z_deriv(vector<Body> o){
            int n = o.size();
            vector<long double> zd(n,0);
            for (int i=0;i<n;i++){
                zd[i]=o[i].momentum[2]/o[i].mass;
            }
            return zd;
        }
        vector<Body> k1(vector<Body>& o,long double h){
            int n = o.size();
            vector<long double> xd1 = x_deriv(o);
            vector<long double> yd1 = y_deriv(o);
            vector<long double> zd1 = z_deriv(o);
            vector<long double> pxd1 = px_deriv(o);
            vector<long double> pyd1 = py_deriv(o);
            vector<long double> pzd1 = pz_deriv(o);

            vector<Body> k1_body = o;
            for (int i=0;i<n;i++){
                k1_body[i].position[0] += h*xd1[i]/2; //x2
                k1_body[i].position[1] += h*yd1[i]/2; //y2
                k1_body[i].position[2] += h*zd1[i]/2; //z2
                k1_body[i].momentum[0] += h*pxd1[i]/2; //px2
                k1_body[i].momentum[1] += h*pyd1[i]/2; //py2
                k1_body[i].momentum[2] += h*pzd1[i]/2; //py2
            }
            return k1_body;
        }
        vector<Body> k2(vector<Body>& o, long double h, vector<Body> k1_body){
            int n = o.size();
            vector<long double> xd2 = x_deriv(k1_body);
            vector<long double> yd2 = y_deriv(k1_body);
            vector<long double> zd2 = z_deriv(k1_body);
            vector<long double> pxd2 = px_deriv(k1_body);
            vector<long double> pyd2 = py_deriv(k1_body);
            vector<long double> pzd2 = pz_deriv(k1_body);
            vector<Body> k2_body = o;
            for (int i=0;i<n;i++){
                k2_body[i].position[0] += h*xd2[i]/2; //x3
                k2_body[i].position[1] += h*yd2[i]/2; //y3
                k2_body[i].position[2] += h*zd2[i]/2; //z3
                k2_body[i].momentum[0] += h*pxd2[i]/2; //px3
                k2_body[i].momentum[1] += h*pyd2[i]/2; //py3
                k2_body[i].momentum[2] += h*pzd2[i]/2; //pz3
            }
            return k2_body;
        }
        vector<Body> k3(vector<Body>& o, long double h, vector<Body> k2_body){
            int n = o.size();
            vector<long double> xd3 = x_deriv(k2_body);
            vector<long double> yd3 = y_deriv(k2_body);
            vector<long double> zd3 = z_deriv(k2_body);
            vector<long double> pxd3 = px_deriv(k2_body);
            vector<long double> pyd3 = py_deriv(k2_body);
            vector<long double> pzd3 = pz_deriv(k2_body);
            vector<Body> k3_body = o;
            for (int i=0;i<n;i++){
                k3_body[i].position[0] += h*xd3[i]; //x4
                k3_body[i].position[1] += h*yd3[i]; //y4
                k3_body[i].position[2] += h*zd3[i]; //z4
                k3_body[i].momentum[0] += h*pxd3[i]; //px4
                k3_body[i].momentum[1] += h*pyd3[i]; //py4
                k3_body[i].momentum[2] += h*pzd3[i]; //pz4
            }
            return k3_body;
        }
        vector<Body> k4(vector<Body>& o, long double h, vector<Body> k3_body){
            int n = o.size();
            vector<long double> xd4 = x_deriv(k3_body);
            vector<long double> yd4 = y_deriv(k3_body);
            vector<long double> zd4 = z_deriv(k3_body);
            vector<long double> pxd4 = px_deriv(k3_body);
            vector<long double> pyd4 = py_deriv(k3_body);
            vector<long double> pzd4 = pz_deriv(k3_body);
            vector<Body> k4_body = o;
            for (int i=0;i<n;i++){
                k4_body[i].position[0] += h*xd4[i]; 
                k4_body[i].position[1] += h*yd4[i]; 
                k4_body[i].position[2] += h*zd4[i];
                k4_body[i].momentum[0] += h*pxd4[i];
                k4_body[i].momentum[1] += h*pyd4[i];
                k4_body[i].momentum[2] += h*pzd4[i];
            }
            return k4_body;
        }
    public:
        void rk4(vector<Body>& o, long double h){
            int n = o.size();
            vector<Body> k1_body = k1(o,h);
            vector<Body> k2_body = k2(o,h,k1_body);
            vector<Body> k3_body = k3(o,h,k2_body);
            vector<Body> k4_body = k4(o,h,k3_body);
            for (int i=0;i<n;i++){
                o[i].position[0] += (2*k1_body[i].position[0]+4*k2_body[i].position[0]+2*k3_body[i].position[0]+k4_body[i].position[0]-9*o[i].position[0])/6;
                o[i].position[1] += (2*k1_body[i].position[1]+4*k2_body[i].position[1]+2*k3_body[i].position[1]+k4_body[i].position[1]-9*o[i].position[1])/6;
                o[i].position[2] += (2*k1_body[i].position[2]+4*k2_body[i].position[2]+2*k3_body[i].position[2]+k4_body[i].position[2]-9*o[i].position[2])/6;
                o[i].momentum[0] += (2*k1_body[i].momentum[0]+4*k2_body[i].momentum[0]+2*k3_body[i].momentum[0]+k4_body[i].momentum[0]-9*o[i].momentum[0])/6;
                o[i].momentum[1] += (2*k1_body[i].momentum[1]+4*k2_body[i].momentum[1]+2*k3_body[i].momentum[1]+k4_body[i].momentum[1]-9*o[i].momentum[1])/6;
                o[i].momentum[2] += (2*k1_body[i].momentum[2]+4*k2_body[i].momentum[2]+2*k3_body[i].momentum[2]+k4_body[i].momentum[2]-9*o[i].momentum[2])/6; 
            }
            collision(o,0.4);
        }

        long double totenergy(vector<Body> o){
            int n = o.size();
            long double ke = 0;
            long double pe =0;
            for (int i=0;i<n;i++){
                ke += (pow(o[i].momentum[0],2)+pow(o[i].momentum[1],2)+pow(o[i].momentum[2],2))/(2*o[i].mass);
            }
            for (int i=0;i<n-1;i++){
                for (int j=i+1;j<n;j++){
                    long double xr = o[j].position[0] - o[i].position[0];
                    long double yr = o[j].position[1] - o[i].position[1];
                    long double zr = o[j].position[2] - o[i].position[2];
                    pe -= (G*o[i].mass*o[j].mass)/sqrtl(xr*xr+yr*yr+zr*zr);
                }
            }
            long double totale = ke + pe;
            return totale;
        }

        vector<long double> totmomentum(vector<Body> o){
            int n = o.size();
            long double px =0;
            long double py =0;
            long double pz =0;
            for (int i=0;i<n;i++){
                px += o[i].momentum[0];
                py += o[i].momentum[1];
                pz += o[i].momentum[2];
            }
            return {px,py,pz};
        }
};

int main(){
    vector<Body> o ={Body(1988500e24,696000e3,-1.239974922284681E+09,-3.575289179956036E+08,3.187652087549341E+07,7.355252437281626,-1.330967142719719E1,-5.509223105502978E-02),
    Body(3.302e23,2440e3,-5.468476275907741E+010,-3.797389382366690E+010,1.860013342376858E+09,1.805937653404739E+04,-3.769137328626643E+04,-4.734952989284444E+03),
    Body(48.685e23,6051.84e3,3.091837907949718E+010,1.025646766503651E+11,4.102134207988381E+08,-3.353555387985905E+04,1.026385491922525E+04,2.076534282628684E+03),
    Body(5.97219e24,6371.01e3,1.289707090410556E+11,7.180956016858132E+10,2.742526014676318E+07,-1.492708726900149E+04,2.592116171618013E+04,-1.741940167978129),
    Body(6.4171e23,3389.92e3,-1.754080483637147E+011,-1.579333162024269E+11,1.001739413205221E+09,1.717113714481583E+03,-1.590955897727394E+04,-7.542166771319083E+02)};
    long double t_start = 0;
    long double t_end = 3600*24*365*10;
    long double h = 3600;
    solver gravity_solver;
    ofstream file1; file1.open("foo2.csv");
    for(long double t=t_start;t<=t_end;t+=h){
        gravity_solver.rk4(o,h);
        for (int j =0;j<o.size();j++){
            file1<<o[j].position[0]<<","<<o[j].position[1]<<","<<o[j].position[2]<<",";
        }
        file1<<endl;
    }   
    file1.close();
    cout<<"Completed"<<endl;
    return 0;
}