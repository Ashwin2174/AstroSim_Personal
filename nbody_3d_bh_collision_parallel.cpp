#include <bits/stdc++.h>
#include <vector>
#include <omp.h>
using namespace std;
const long double G = 1;

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

class Node{
            public:
                long double mass = 0; // total mass of node
                vector<long double> com; // centre of mass of node
                vector<Body*> particles; //ptrs to particles inside the node
                vector<long double> centre; // centre of node
                long double s; // side length
                vector<Node*> child; // ptrs to child nodes
                Node(long double x, long double y, long double z, long double side){
                    centre = {x,y,z};
                    s = side;
                }
};

class BarnesHut {
    private:
        void Divide(Node*& node){
            long double x = node->centre[0];
            long double y = node->centre[1];
            long double z = node->centre[2];
            long double s = node->s;
            Node* btr = new Node(x+s/4,y+s/4,z-s/4,s/2); // b-bottom t-top l-left r-right
            Node* btl = new Node(x-s/4,y+s/4,z-s/4,s/2);
            Node* bbl = new Node(x-s/4,y-s/4,z-s/4,s/2);
            Node* bbr = new Node(x+s/4,y-s/4,z-s/4,s/2);
            Node* ttr = new Node(x+s/4,y+s/4,z+s/4,s/2);
            Node* ttl = new Node(x-s/4,y+s/4,z+s/4,s/2);
            Node* tbl = new Node(x-s/4,y-s/4,z+s/4,s/2);
            Node* tbr = new Node(x+s/4,y-s/4,z+s/4,s/2);
            node->child = {btr,btl,bbl,bbr,ttr,ttl,tbl,tbr};
        }

        int whichOctant(Node*& node, Body& body){
            long double xr = body.position[0] - node->centre[0];
            long double yr = body.position[1] - node->centre[1];
            long double zr = body.position[2] - node->centre[2];
            if (xr>0 && yr>0 && zr<0){return 0;}
            else if (xr<0 && yr>0 && zr<0){return 1;}
            else if (xr<0 && yr<0 && zr<0){return 2;}
            else if (xr>0 && yr<0 && zr<0){return 3;}
            else if (xr>0 && yr>0 && zr>0){return 4;}
            else if (xr<0 && yr>0 && zr>0){return 5;}
            else if (xr<0 && yr<0 && zr>0){return 6;}
            else {return 7;}
        }

        void Octinsert(Node*& node, Body& body){
            if (CurrentDepth < MaxRecursionDepth){
                if (node->particles.size()>1){
                    node->com[0] = (node->mass*node->com[0] + body.mass*body.position[0])/(node->mass+body.mass);
                    node->com[1] = (node->mass*node->com[1] + body.mass*body.position[1])/(node->mass+body.mass);
                    node->com[2] = (node->mass*node->com[2] + body.mass*body.position[2])/(node->mass+body.mass);
                    node->mass += body.mass;
                    node->particles.push_back(&body);
                    CurrentDepth += 1;
                    Octinsert(node->child[whichOctant(node,body)],body);
                }
                else if (node->particles.size()==1){
                    Divide(node);
                    CurrentDepth += 1;
                    Octinsert(node->child[whichOctant(node,*(node->particles[0]))],*(node->particles[0]));
                    node->com[0] = (node->mass*node->com[0] + body.mass*body.position[0])/(node->mass+body.mass);
                    node->com[1] = (node->mass*node->com[1] + body.mass*body.position[1])/(node->mass+body.mass);
                    node->com[2] = (node->mass*node->com[2] + body.mass*body.position[2])/(node->mass+body.mass);
                    node->mass += body.mass;
                    node->particles.push_back(&body);
                    Octinsert(node->child[whichOctant(node,body)],body);
                }
                else {
                    node->com = body.position;
                    node->mass = body.mass;
                    node->particles = {&body};
                }
            }
            else{
                if (node->particles.size()>=1){
                    node->com[0] = (node->mass*node->com[0] + body.mass*body.position[0])/(node->mass+body.mass);
                    node->com[1] = (node->mass*node->com[1] + body.mass*body.position[1])/(node->mass+body.mass);
                    node->com[2] = (node->mass*node->com[2] + body.mass*body.position[2])/(node->mass+body.mass);
                    node->mass += body.mass;
                    node->particles.push_back(&body);
                }
                else{
                    node->com = body.position;
                    node->mass = body.mass;
                    node->particles = {&body};
                }
            }
        }

        void DeleteUnusedNodes(Node*& root){
            for (auto it=root->child.begin();it!=root->child.end();){
                if ((*it)->particles.size()==0){
                    delete *it;
                    *it = NULL;
                    it = root->child.erase(it);
                }
                else if ((*it)->child.size()==0){
                    it++;
                }
                else {
                    DeleteUnusedNodes(*it);
                    it++;
                }
            }
        }

        void OcttreeBuild(Node*& root, vector<Body>& o){
            for (int i=0;i<o.size();i++){
                Octinsert(root,o[i]);
                CurrentDepth = 1;
            }
            DeleteUnusedNodes(root);
        }

        void DeleteOcttree(Node*& root){
            for (auto it=root->child.begin();it!=root->child.end();it++){
                if ((*it)->child.size()==0){
                    delete *it;
                    *it = NULL;
                }
                else {
                    DeleteOcttree(*it);
                }
            }
            delete root;
            root = NULL;
        }

        vector<long double> TreeForce(Node*& node, Body& body){
            long double fx=0, fy=0, fz=0;
            long double r = sqrt(pow((node->com[0]-body.position[0]),2)+pow((node->com[1]-body.position[1]),2)+pow((node->com[2]-body.position[2]),2));
            long double theta = 0.5;
            if (node->s < theta*r){
                fx = G*body.mass*node->mass*(node->com[0]-body.position[0])/pow(r,3);
                fy = G*body.mass*node->mass*(node->com[1]-body.position[1])/pow(r,3);
                fz = G*body.mass*node->mass*(node->com[2]-body.position[2])/pow(r,3);
            }
            else{
                if (node->child.size()>0){
                    for (auto it=node->child.begin();it!=node->child.end();it++){
                        vector<long double> ChildForce = TreeForce(*it,body);
                        fx+=ChildForce[0];
                        fy+=ChildForce[1];
                        fz+=ChildForce[2];
                    } 
                }
                else{
                    for (int i=0;i<node->particles.size();i++){
                        if (node->particles[i] != &body){
                            long double r2 = pow((node->particles[i]->position[0]-body.position[0]),2)+pow((node->particles[i]->position[1]-body.position[1]),2)+pow((node->particles[i]->position[2]-body.position[2]),2);
                            fx += G*body.mass*node->particles[i]->mass*(node->particles[i]->position[0]-body.position[0])/(pow(r2,1.5));
                            fy += G*body.mass*node->particles[i]->mass*(node->particles[i]->position[1]-body.position[1])/(pow(r2,1.5));
                            fz += G*body.mass*node->particles[i]->mass*(node->particles[i]->position[2]-body.position[2])/(pow(r2,1.5));  
                        }
                    }
                }
            }
            return {fx,fy,fz};
        }

        bool intersect(vector<long double>& cube, vector<long double>& sphere, long double r, long double side){
            vector<long double> C1 = {cube[0]-side/2,cube[1]-side/2,cube[2]-side/2};
            vector<long double> C2 = {cube[0]+side/2,cube[1]+side/2,cube[2]+side/2};
            long double dist_squared = r * r;
            if (sphere[0] < C1[0]){
                dist_squared -= (sphere[0] - C1[0])*(sphere[0] - C1[0]);
            } 
            else if (sphere[0] > C2[0]){
                dist_squared -= (sphere[0] - C2[0])*(sphere[0] - C2[0]);
            } 
            if (sphere[1] < C1[1]){
                dist_squared -= (sphere[1] - C1[1])*(sphere[1] - C1[1]);
            } 
            else if (sphere[1] > C2[1]){
                dist_squared -= (sphere[1] - C2[1])*(sphere[1] - C2[1]);
            } 
            if (sphere[2] < C1[2]){
                dist_squared -= (sphere[2] - C1[2])*(sphere[2] - C1[2]);
            } 
            else if (sphere[2] > C2[2]){
                dist_squared -= (sphere[2] - C2[2])*(sphere[2] - C2[2]);
            } 
            return dist_squared > 0;
        }

        void query(Node*& node, Body& body, long double maxradius, vector<Body*>& neighbours){
            if (intersect(node->centre,body.position,maxradius+body.radius,node->s)){
                if (node->child.size()>0){
                    for (auto it=node->child.begin();it!=node->child.end();it++){
                        query(*it,body,maxradius,neighbours);
                    }
                }
                else{
                    for (int i =0;i<node->particles.size();i++){
                        long double d2 = pow((node->particles[i]->position[0]-body.position[0]),2)+pow((node->particles[i]->position[1]-body.position[1]),2)+pow((node->particles[i]->position[2]-body.position[2]),2);
                        if (d2 <= pow(body.radius+node->particles[i]->radius,2)){
                            neighbours.push_back(node->particles[i]);
                        }
                    }
                }
            }
        }


    public:
        int CurrentDepth = 1;
        int MaxRecursionDepth = 100;
        long double CoefficientOfRestitution = 0.5;
        long double MaxRadius = 0;

        vector<long double> RootCenter;
        long double RootSide;

        long double ComputeRootSide(vector<Body>& o){
            long double side = 0;
            for (int i=0;i<o.size();i++){
                side = max({side,2*abs(o[i].position[0]),2*abs(o[i].position[1]),2*abs(o[i].position[2])});
            }
            return side;
        }

        long double MaximumRadius(vector<Body>& o){
            long double mxradius = 0;
            for (int i=0;i<o.size();i++){
                mxradius = max(o[i].radius,mxradius);
            }
            return mxradius;
        }

        BarnesHut(vector<long double> center){
            RootCenter = center;
        }

        vector<vector<long double>> BarnesHutForces(vector<Body>& o){
            RootSide = ComputeRootSide(o);
            Node* root = new Node(RootCenter[0],RootCenter[1],RootCenter[2],RootSide);
            OcttreeBuild(root,o);
            vector<long double> fx(o.size(),0), fy(o.size(),0), fz(o.size(),0);
            #pragma omp parallel for
            for (int i=0;i<o.size();i++){
                vector<long double> ForceOnBody = TreeForce(root,o[i]);
                fx[i] = ForceOnBody[0];
                fy[i] = ForceOnBody[1];
                fz[i] = ForceOnBody[2];
            }
            DeleteOcttree(root);
            return {fx,fy,fz};
        }

        void OcttreeCollision(vector<Body>& o){
            RootSide = ComputeRootSide(o);
            Node* root = new Node(RootCenter[0],RootCenter[1],RootCenter[2],RootSide);
            OcttreeBuild(root,o);
            int CollisionCount = 0;
            for (int i =0;i<o.size();i++){
                vector<Body*> neighbours;
                query(root,o[i],MaxRadius,neighbours);
                for (int j=0;j<neighbours.size();j++){
                    if (neighbours[j] != &o[i]){
                        CollisionCount+=1;
                        long double delx = neighbours[j]->position[0] - o[i].position[0];
                        long double dely = neighbours[j]->position[1] - o[i].position[1];
                        long double delz = neighbours[j]->position[2] - o[i].position[2];
                        long double r = sqrt(delx*delx+dely*dely+delz*delz);
                        long double k = o[i].radius+neighbours[j]->radius-r;
                        long double nx = delx/r;
                        long double ny = dely/r;
                        long double nz = delz/r;
                        long double v_impact = nx*(o[i].momentum[0]/o[i].mass-neighbours[j]->momentum[0]/neighbours[j]->mass) + ny*(o[i].momentum[1]/o[i].mass-neighbours[j]->momentum[1]/neighbours[j]->mass) + nz*(o[i].momentum[2]/o[i].mass-neighbours[j]->momentum[2]/neighbours[j]->mass);
                        long double reduced_mass = (o[i].mass*neighbours[j]->mass)/(o[i].mass+neighbours[j]->mass);
                        long double impulse = (1+CoefficientOfRestitution)*reduced_mass*v_impact;
                        o[i].momentum[0] -= impulse*nx;
                        o[i].momentum[1] -= impulse*ny;
                        o[i].momentum[2] -= impulse*nz;
                        neighbours[j]->momentum[0] += impulse*nx;
                        neighbours[j]->momentum[1] += impulse*ny;
                        neighbours[j]->momentum[2] += impulse*nz;
                        o[i].position[0] -= 0.51*k*nx;
                        neighbours[j]->position[0] += 0.51*k*nx;
                        neighbours[j]->position[1] += 0.51*k*ny;
                        o[i].position[1] -= 0.51*k*ny;
                        o[i].position[2] -= 0.51*k*nz;
                        neighbours[j]->position[2] += 0.51*k*nz; 
                    }
                }
            }
            if (CollisionCount>0){
                std::cout<<"Collisions : "<<CollisionCount<<endl;
            }
            DeleteOcttree(root);
        }

};


class solver: public BarnesHut {
    private:
        vector<long double> x_deriv(vector<Body> o){
            int n = o.size();
            vector<long double> xd(n,0);
            #pragma omp parallel for
            for (int i=0;i<n;i++){
                xd[i]=o[i].momentum[0]/o[i].mass;
            }
            return xd;
        }
        vector<long double> y_deriv(vector<Body> o){
            int n = o.size();
            vector<long double> yd(n,0);
            #pragma omp parallel for
            for (int i=0;i<n;i++){
                yd[i]=o[i].momentum[1]/o[i].mass;
            }
            return yd;
        }
        vector<long double> z_deriv(vector<Body> o){
            int n = o.size();
            vector<long double> zd(n,0);
            #pragma omp parallel for
            for (int i=0;i<n;i++){
                zd[i]=o[i].momentum[2]/o[i].mass;
            }
            return zd;
        }
        vector<Body> k1(vector<Body>& o,long double h){
            vector<vector<long double>> forces = BarnesHutForces(o);
            int n = o.size();
            vector<long double> xd1 = x_deriv(o);
            vector<long double> yd1 = y_deriv(o);
            vector<long double> zd1 = z_deriv(o);
            vector<long double> pxd1 = forces[0];
            vector<long double> pyd1 = forces[1];
            vector<long double> pzd1 = forces[2];

            vector<Body> k1_body = o;
            #pragma omp parallel for
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
        vector<Body> k2(vector<Body>& o, long double h, vector<Body>& k1_body){
            vector<vector<long double>> forces = BarnesHutForces(k1_body);
            int n = o.size();
            vector<long double> xd2 = x_deriv(k1_body);
            vector<long double> yd2 = y_deriv(k1_body);
            vector<long double> zd2 = z_deriv(k1_body);
            vector<long double> pxd2 = forces[0];
            vector<long double> pyd2 = forces[1];
            vector<long double> pzd2 = forces[2];
            vector<Body> k2_body = o;
            #pragma omp parallel for
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
        vector<Body> k3(vector<Body>& o, long double h, vector<Body>& k2_body){
            vector<vector<long double>> forces = BarnesHutForces(k2_body);
            int n = o.size();
            vector<long double> xd3 = x_deriv(k2_body);
            vector<long double> yd3 = y_deriv(k2_body);
            vector<long double> zd3 = z_deriv(k2_body);
            vector<long double> pxd3 = forces[0];
            vector<long double> pyd3 = forces[1];
            vector<long double> pzd3 = forces[2];
            vector<Body> k3_body = o;
            #pragma omp parallel for
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
        vector<Body> k4(vector<Body>& o, long double h, vector<Body>& k3_body){
            vector<vector<long double>> forces = BarnesHutForces(k3_body);
            int n = o.size();
            vector<long double> xd4 = x_deriv(k3_body);
            vector<long double> yd4 = y_deriv(k3_body);
            vector<long double> zd4 = z_deriv(k3_body);
            vector<long double> pxd4 = forces[0];
            vector<long double> pyd4 = forces[1];
            vector<long double> pzd4 = forces[2];
            vector<Body> k4_body = o;
            #pragma omp parallel for
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
        solver(vector<long double> center):BarnesHut(center){}
        void rk4(vector<Body>& o, long double h){
            int n = o.size();
            vector<Body> k1_body = k1(o,h);
            vector<Body> k2_body = k2(o,h,k1_body);
            vector<Body> k3_body = k3(o,h,k2_body);
            vector<Body> k4_body = k4(o,h,k3_body);
            #pragma omp parallel for
            for (int i=0;i<n;i++){
                o[i].position[0] += (2*k1_body[i].position[0]+4*k2_body[i].position[0]+2*k3_body[i].position[0]+k4_body[i].position[0]-9*o[i].position[0])/6;
                o[i].position[1] += (2*k1_body[i].position[1]+4*k2_body[i].position[1]+2*k3_body[i].position[1]+k4_body[i].position[1]-9*o[i].position[1])/6;
                o[i].position[2] += (2*k1_body[i].position[2]+4*k2_body[i].position[2]+2*k3_body[i].position[2]+k4_body[i].position[2]-9*o[i].position[2])/6;
                o[i].momentum[0] += (2*k1_body[i].momentum[0]+4*k2_body[i].momentum[0]+2*k3_body[i].momentum[0]+k4_body[i].momentum[0]-9*o[i].momentum[0])/6;
                o[i].momentum[1] += (2*k1_body[i].momentum[1]+4*k2_body[i].momentum[1]+2*k3_body[i].momentum[1]+k4_body[i].momentum[1]-9*o[i].momentum[1])/6;
                o[i].momentum[2] += (2*k1_body[i].momentum[2]+4*k2_body[i].momentum[2]+2*k3_body[i].momentum[2]+k4_body[i].momentum[2]-9*o[i].momentum[2])/6; 
            }
        }

        void leapfrog(vector<Body>& o, long double h){
            int n = o.size();
            for (int i=0;i<n;i++){
                o[i].position[0] += 0.5*h*o[i].momentum[0]/o[i].mass; 
                o[i].position[1] += 0.5*h*o[i].momentum[1]/o[i].mass; 
                o[i].position[2] += 0.5*h*o[i].momentum[2]/o[i].mass;
            }
            vector<vector<long double>> forces = BarnesHutForces(o);
            vector<long double> pxd3 = forces[0];
            vector<long double> pyd3 = forces[1];
            vector<long double> pzd3 = forces[2];
            for (int i=0;i<n;i++){
                o[i].momentum[0] += h*pxd3[i]; //px4
                o[i].momentum[1] += h*pyd3[i]; //py4
                o[i].momentum[2] += h*pzd3[i]; //pz4
                o[i].position[0] += 0.5*h*o[i].momentum[0]/o[i].mass; 
                o[i].position[1] += 0.5*h*o[i].momentum[1]/o[i].mass; 
                o[i].position[2] += 0.5*h*o[i].momentum[2]/o[i].mass;
            }
        }

        long double totenergy(vector<Body>& o){
            int n = o.size();
            long double ke = 0;
            long double pe =0;
            #pragma omp parallel for reduction (+:ke)
            for (int i=0;i<n;i++){
                ke += (pow(o[i].momentum[0],2)+pow(o[i].momentum[1],2)+pow(o[i].momentum[2],2))/(2*o[i].mass);
            }
            #pragma omp parallel for reduction (+:pe)
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

        vector<long double> totmomentum(vector<Body>& o){
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
    long double t_start = 0;
    long double t_end = 1;
    long double h = 0.04;

    solver gravity_solver = solver({0,0,0});
    ifstream input("galaxy.tab");
    vector<Body> o;
    for (int i =0;i<81920;i++){
        long double m,x,y,z,vx,vy,vz,id,r;
        input>>m>>x>>y>>z>>vx>>vy>>vz; 
        if (i%1==0){
        o.push_back(Body(m,0,x,y,z,vx,vy,vz));
        }
    }
    cout<<o.size()<<endl;
    cout<<"-------------------------------------------------------------------"<<endl;
    ofstream output; output.open("galaxsim.csv");

    //gravity_solver.MaxRadius = gravity_solver.MaximumRadius(o);
    for(long double t=t_start;t<=t_end;t+=h){
        auto it_start = std::chrono::high_resolution_clock::now();
        gravity_solver.rk4(o,h);
        //gravity_solver.OcttreeCollision(o);
        if (int(t/h)%1==0){
            for (int j =0;j<o.size();j++){
                if (j!=o.size()-1){
                    output<<o[j].position[0]<<","<<o[j].position[1]<<","<<o[j].position[2]<<",";
                }
                else {
                    output<<o[j].position[0]<<","<<o[j].position[1]<<","<<o[j].position[2]<<'\n';
                } 
            }
        }
        auto it_end = std::chrono::high_resolution_clock::now();
        auto it_time = std::chrono::duration_cast<std::chrono::milliseconds>(it_end - it_start);
        cout<<"Iteration "<<int(t/h)<<" : "<<it_time.count()<<" ms"<<endl;
        //vector<long double> tm = gravity_solver.totmomentum(o);
        //cout<<"Energy : "<<gravity_solver.totenergy(o)<<" J Momentum : "<<tm[0]<<" "<<tm[1]<<" "<<tm[2]<<" kgm/s"<<endl;
        cout<<"-------------------------------------------------------------------"<<endl;
    }   
    output.close();
    std::cout<<"Completed"<<endl;
    return 0;
}