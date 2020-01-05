#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <array>

namespace patch
{
/* patch to make to_string work */
template < typename T > std::string to_string( const T& n )
{
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
}
}

using namespace std;
struct coord
{
    long double x;
    long double y;
    long double z;
    bool operator==(const coord& a) const
    {
          return (x == a.x && y == a.y && z == a.z);
    }
};

struct medians
{
    coord left;
    coord right;
    coord bottom;
    coord top;
    coord back;
    coord front;
};

struct vertices
{
    coord a;
    coord b;
    coord c;
    coord d;
    coord e;
    coord f;
    coord g;
    coord h;
};

struct faces
{
    std::string left,right,bottom,top,back,front;
};



class read_compute
{
public:
    int tot_grids;
    vector <int> imax, jmax, kmax;
    std::vector <std::string> grid_names;
    std::vector <vertices> all_vertices;
    std::vector <faces> all_faces;
    std::vector<std::array<int, 6>> mpi_class;
    void read_grid_names(std::ostream& block);//read file_names to get grids
    void compute_grid_medians(std::string str, std::ostream& block );
    void compute_interfaces();
    void write_to_file();
    void write_to_mapfile(std::ostream& map, int b1,int f1,int s11,int s12,int e11,
    int e12,int b2,int f2,int s21,int s22,int e21,int e22, bool dir_swap);
    void check_interface(int b1, int b2, int fb, int fc, 
        coord v11, coord v12, coord v13, coord v14,
        coord v21, coord v22, coord v23, coord v24, std::ostream& map);
    std::string layout_comment(std::string str);
};

void read_compute::write_to_file()
{

    string str = "layout.md",comnt_str;
    char *c = &str[0u];
    std::fstream output_file (c,fstream::out);
    if(!output_file.is_open())
    {
        cout << "Error: Couldn't open output file" << endl;
        return;
    }
    comnt_str = "BLOCK LAYOUT FILE";
    output_file << layout_comment(comnt_str);
    comnt_str = "==========================";
    output_file << layout_comment(comnt_str);
    comnt_str = "NUMBER OF PROCESSES";
    output_file << layout_comment(comnt_str);
    output_file << tot_grids << "\n";
    comnt_str = "NUMBER OF ENTRIES PER PROCESS";
    output_file << layout_comment(comnt_str);
    output_file << "9\n";
    comnt_str = "PROCESS_NO GRID BC_FILE IMIN IMAX JMIN JMAX KMIN KMAX";
    output_file << layout_comment(comnt_str);
    comnt_str = "===================================";
    output_file << layout_comment(comnt_str);
    for (int i =0; i< tot_grids; i++)
    {
    comnt_str = "PROCESS "+ patch::to_string(i);
    output_file << layout_comment(comnt_str);
    output_file <<setfill('0') <<setw(2) <<i << "  " << grid_names[i].substr(13,-1)<<"  "<<"bc_" << setw(2) << setfill('0') << i<<".md" <<"  ";
    output_file << setw(4) << all_faces[i].left << "  " << setw(4) << all_faces[i].right << "  "<< setw(4) << all_faces[i].bottom << "  ";
    output_file << setw(4) << all_faces[i].top << "  "<< setw(4) << all_faces[i].front << "  "<< setw(4) << all_faces[i].back << "  ";
    output_file << "\n";
    }
    output_file.close();

}

void read_compute::compute_interfaces()
{

    bool dir_swap=false;
    vertices buf1,buf2;
    coord b_Va,b_Vb,b_Vc,b_Vd;
    coord b_Ve,b_Vf,b_Vg,b_Vh;
    coord c_Va,c_Vb,c_Vc,c_Vd;
    coord c_Ve,c_Vf,c_Vg,c_Vh;
    std::ofstream map;
    map.open("mapping.txt");
    map << "#   B1    F1    S1    E1    S2    E2    B2    F2    S1    E1    S2    E2    dir_swap" << endl;
    for(int i =0; i< tot_grids ; i++)
    {
        /*
        imin face a -> c -> e -> g
        imax face b -> d -> f -> h
        jmin face a -> b -> e -> f
        jmax face c -> d -> g -> h
        kmin face a -> b -> c -> d
        kmax face e -> f -> g -> h
        */
        // left face check
//        if(all_faces[i].left == "USER")
//        {
           //cout << "enterd 1";
            buf1 = all_vertices[i];
            for(int j = 0; j< tot_grids; j++)
            {
              if(i!=j){
                buf2 = all_vertices[j];
                b_Va = buf1.a;
                b_Vb = buf1.b;
                b_Vc = buf1.c;
                b_Vd = buf1.d;
                b_Ve = buf1.e;
                b_Vf = buf1.f;
                b_Vg = buf1.g;
                b_Vh = buf1.h;
                c_Va = buf2.a;
                c_Vb = buf2.b;
                c_Vc = buf2.c;
                c_Vd = buf2.d;
                c_Ve = buf2.e;
                c_Vf = buf2.f;
                c_Vg = buf2.g;
                c_Vh = buf2.h;
                //cout << c_buf1.x << " " << c_buf2.x << endl;
                check_interface(i,j,1,1,b_Va,b_Vc,b_Vg,b_Ve,c_Va,c_Vc,c_Vg,c_Ve,map);
                check_interface(i,j,1,2,b_Va,b_Vc,b_Vg,b_Ve,c_Vb,c_Vd,c_Vh,c_Vf,map);
                check_interface(i,j,1,3,b_Va,b_Vc,b_Vg,b_Ve,c_Va,c_Vb,c_Vf,c_Ve,map);
                check_interface(i,j,1,4,b_Va,b_Vc,b_Vg,b_Ve,c_Vc,c_Vd,c_Vh,c_Vg,map);
                check_interface(i,j,1,5,b_Va,b_Vc,b_Vg,b_Ve,c_Va,c_Vb,c_Vd,c_Vc,map);
                check_interface(i,j,1,6,b_Va,b_Vc,b_Vg,b_Ve,c_Ve,c_Vf,c_Vh,c_Vg,map);

                check_interface(i,j,2,1,b_Vb,b_Vd,b_Vh,b_Vf,c_Va,c_Vc,c_Vg,c_Ve,map);
                check_interface(i,j,2,2,b_Vb,b_Vd,b_Vh,b_Vf,c_Vb,c_Vd,c_Vh,c_Vf,map);
                check_interface(i,j,2,3,b_Vb,b_Vd,b_Vh,b_Vf,c_Va,c_Vb,c_Vf,c_Ve,map);
                check_interface(i,j,2,4,b_Vb,b_Vd,b_Vh,b_Vf,c_Vc,c_Vd,c_Vh,c_Vg,map);
                check_interface(i,j,2,5,b_Vb,b_Vd,b_Vh,b_Vf,c_Va,c_Vb,c_Vd,c_Vc,map);
                check_interface(i,j,2,6,b_Vb,b_Vd,b_Vh,b_Vf,c_Ve,c_Vf,c_Vh,c_Vg,map);

                check_interface(i,j,3,1,b_Va,b_Vb,b_Vf,b_Ve,c_Va,c_Vc,c_Vg,c_Ve,map);
                check_interface(i,j,3,2,b_Va,b_Vb,b_Vf,b_Ve,c_Vb,c_Vd,c_Vh,c_Vf,map);
                check_interface(i,j,3,3,b_Va,b_Vb,b_Vf,b_Ve,c_Va,c_Vb,c_Vf,c_Ve,map);
                check_interface(i,j,3,4,b_Va,b_Vb,b_Vf,b_Ve,c_Vc,c_Vd,c_Vh,c_Vg,map);
                check_interface(i,j,3,5,b_Va,b_Vb,b_Vf,b_Ve,c_Va,c_Vb,c_Vd,c_Vc,map);
                check_interface(i,j,3,6,b_Va,b_Vb,b_Vf,b_Ve,c_Ve,c_Vf,c_Vh,c_Vg,map);
                
                check_interface(i,j,4,1,b_Vc,b_Vd,b_Vh,b_Vg,c_Va,c_Vc,c_Vg,c_Ve,map);
                check_interface(i,j,4,2,b_Vc,b_Vd,b_Vh,b_Vg,c_Vb,c_Vd,c_Vh,c_Vf,map);
                check_interface(i,j,4,3,b_Vc,b_Vd,b_Vh,b_Vg,c_Va,c_Vb,c_Vf,c_Ve,map);
                check_interface(i,j,4,4,b_Vc,b_Vd,b_Vh,b_Vg,c_Vc,c_Vd,c_Vh,c_Vg,map);
                check_interface(i,j,4,5,b_Vc,b_Vd,b_Vh,b_Vg,c_Va,c_Vb,c_Vd,c_Vc,map);
                check_interface(i,j,4,6,b_Vc,b_Vd,b_Vh,b_Vg,c_Ve,c_Vf,c_Vh,c_Vg,map);

                check_interface(i,j,5,1,b_Va,b_Vb,b_Vd,b_Vc,c_Va,c_Vc,c_Vg,c_Ve,map);
                check_interface(i,j,5,2,b_Va,b_Vb,b_Vd,b_Vc,c_Vb,c_Vd,c_Vh,c_Vf,map);
                check_interface(i,j,5,3,b_Va,b_Vb,b_Vd,b_Vc,c_Va,c_Vb,c_Vf,c_Ve,map);
                check_interface(i,j,5,4,b_Va,b_Vb,b_Vd,b_Vc,c_Vc,c_Vd,c_Vh,c_Vg,map);
                check_interface(i,j,5,5,b_Va,b_Vb,b_Vd,b_Vc,c_Va,c_Vb,c_Vd,c_Vc,map);
                check_interface(i,j,5,6,b_Va,b_Vb,b_Vd,b_Vc,c_Ve,c_Vf,c_Vh,c_Vg,map);

                check_interface(i,j,6,1,b_Ve,b_Vf,b_Vh,b_Vg,c_Va,c_Vc,c_Vg,c_Ve,map);
                check_interface(i,j,6,2,b_Ve,b_Vf,b_Vh,b_Vg,c_Vb,c_Vd,c_Vh,c_Vf,map);
                check_interface(i,j,6,3,b_Ve,b_Vf,b_Vh,b_Vg,c_Va,c_Vb,c_Vf,c_Ve,map);
                check_interface(i,j,6,4,b_Ve,b_Vf,b_Vh,b_Vg,c_Vc,c_Vd,c_Vh,c_Vg,map);
                check_interface(i,j,6,5,b_Ve,b_Vf,b_Vh,b_Vg,c_Va,c_Vb,c_Vd,c_Vc,map);
                check_interface(i,j,6,6,b_Ve,b_Vf,b_Vh,b_Vg,c_Ve,c_Vf,c_Vh,c_Vg,map);

//                if(b_Va==c_Vb && b_Vc==c_Vd && b_Ve==c_Vf && b_Vg==c_Vh){
//                  all_faces[i].left =patch::to_string(j);
//                  all_faces[j].right =patch::to_string(i);
//                  swap=false;
//                  write_to_mapfile(map,i,1,1,jmax[i],1,kmax[i],j,3,1,jmax[j],1,kmax[j],swap);
//                  break;
//                }
//                if(b_Va==c_Vd && b_Vc==c_Vh && b_Ve==c_Va && b_Vg==c_Vf){
//                  all_faces[i].left =patch::to_string(j);
//                  all_faces[j].right =patch::to_string(i);
//                  swap=true;
//                  write_to_mapfile(map,i,1,1,jmax[i],1,kmax[i],j,3,1,jmax[j],kmax[j],1,swap);
//                  break;
//                }
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//
//                    continue;
//                }
//                c_buf1 = buf1.c;
//                c_buf2 = buf2.d;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                //cout << "enterd 10";
//                    continue;
//                }
//                c_buf1 = buf1.e;
//                c_buf2 = buf2.f;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.g;
//                c_buf2 = buf2.h;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                all_faces[i].left =patch::to_string(j);
//                all_faces[j].right =patch::to_string(i);
//                break;
              }
 //           }
        }

//        // top face check
//        if(all_faces[i].bottom == "USER")
//        {
//            buf1 = all_vertices[i];
//            for(int j = 0; j< tot_grids; j++)
//            {
//                buf2 = all_vertices[j];
//                c_buf1 = buf1.a;
//                c_buf2 = buf2.c;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.b;
//                c_buf2 = buf2.d;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.e;
//                c_buf2 = buf2.g;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.f;
//                c_buf2 = buf2.h;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                all_faces[i].bottom =patch::to_string(j);
//                all_faces[j].top =patch::to_string(i);
//                break;
//            }
//        }
//
//        // front face check
//        if(all_faces[i].front == "USER")
//        {
//            buf1 = all_vertices[i];
//            for(int j = 0; j< tot_grids; j++)
//            {
//                buf2 = all_vertices[j];
//                c_buf1 = buf1.e;
//                c_buf2 = buf2.a;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.g;
//                c_buf2 = buf2.c;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.f;
//                c_buf2 = buf2.b;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                c_buf1 = buf1.h;
//                c_buf2 = buf2.d;
//                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
//                {
//                    continue;
//                }
//                all_faces[i].front =patch::to_string(j);
//                all_faces[j].back =patch::to_string(i);
//                break;
//            }
//        }
//
    }
    map.close();

}
void read_compute::read_grid_names(std::ostream& block)
{
    std::string str_buf;
    int count =0;
    std::ifstream file_buf;
    file_buf.open("grid_names");
    while(file_buf >> str_buf)
    {
        grid_names.push_back(str_buf);
        count++;
    }
    tot_grids = count;
    block << "MBLK 	"<<tot_grids << endl;
}

void read_compute::compute_grid_medians(std::string str, std::ostream& block)
{
    char *cr = &str[0u];
    std::ifstream grid;
    grid.open(cr);
    cout << cr << endl;
    int imx =51,jmx =51, kmx =51;
    coord a,b,c,d,e,f,g,h,buf;
    medians buf_med;
    grid >> imx >> jmx >> kmx;
    imax.push_back(imx);
    jmax.push_back(jmx);
    kmax.push_back(kmx);
    block << imx << "	" << jmx << "	" << kmx << endl;
    grid >> a.x >> a.y >> a.z;
    for(int i = 2; i < imx; i++)
        grid >> buf.x >> buf.y >> buf.z;
    grid >> b.x >> b.y >> b.z;
    for(int j = 2; j < jmx; j++)
    {
        for(int i =1 ; i<=imx ; i++)
        {
            grid >> buf.x >> buf.y >> buf.z;
        }
    }
    grid >> c.x >> c.y >> c.z;
    for(int i = 2; i < imx; i++)
        grid >> buf.x >> buf.y >> buf.z;
    grid >> d.x >> d.y >> d.z;

    //skip in k direction
    for (int k =2; k< kmx; k++)
        for(int j = 1; j <= jmx; j++)
        {
            for(int i =1 ; i<=imx ; i++)
            {
                grid >> buf.x >> buf.y >> buf.z;
            }
        }

    grid >> e.x >> e.y >> e.z;
    for(int i = 2; i < imx; i++)
        grid >> buf.x >> buf.y >> buf.z;
    grid >> f.x >> f.y >> f.z;
    for(int j = 2; j < jmx; j++)
    {
        for(int i =1 ; i<=imx ; i++)
        {
            grid >> buf.x >> buf.y >> buf.z;
        }
    }
    grid >> g.x >> g.y >> g.z;
    //cout << g.x << " " << g.y << " "<< g.z << endl;
    for(int i = 2; i < imx; i++)
        grid >> buf.x >> buf.y >> buf.z;
    grid >> h.x >> h.y >> h.z;

    // store vertices
    vertices buf_vert;
    buf_vert.a = a;
    buf_vert.b = b;
    buf_vert.c = c;
    buf_vert.d = d;
    buf_vert.e = e;
    buf_vert.f = f;
    buf_vert.g = g;
    buf_vert.h = h;
    all_vertices.push_back(buf_vert);
    // store faces
    faces buf_face;
    buf_face.back = "USER";
    buf_face.front = "USER";
    buf_face.left = "USER";
    buf_face.right = "USER";
    buf_face.top = "USER";
    buf_face.bottom = "USER";
    all_faces.push_back(buf_face);
    //j10 : faces types for mpi sequence
    std::array<int, 6> buf_mpi_class={-1,-1,-1,-1,-1,-1};
    mpi_class.push_back(buf_mpi_class);
    grid.close();
}
std::string read_compute::layout_comment(string str)
{
return "## " + str + "\n";
}

void read_compute::write_to_mapfile(std::ostream& map, int b1,int f1,int s11,int s12,int e11,
    int e12,int b2,int f2,int s21,int s22,int e21,int e22, bool dir_swap){
map<<"  " << setw(4) << b1 <<"  "<<setw(4) << f1  <<"  "<<setw(4) << s11<<"  "<<setw(4) << s12 
   <<"  " << setw(4) << e11<<"  "<<setw(4) << e12 <<"  "<<setw(4) << b2 <<"  "<<setw(4) << f2
   <<"  " << setw(4) << s21<<"  "<<setw(4) << s22 <<"  "<<setw(4) << e21<<"  "<<setw(4) << e22
   <<"    "<<dir_swap << "  " << mpi_class[b1][f1]<<endl;
}

void read_compute::check_interface(int i, int j, int f1, int f2, coord v11, coord v12, coord v13, coord v14,
        coord v21, coord v22, coord v23, coord v24, std::ostream& map){
  bool dir_swap=false;
  bool found=false;
  int  type=1;
  int  s11=1,s12=1,s21=1,s22=1;
  int  e11,e12,e21,e22;
  if(v11==v21 && v12==v22 && v13==v23 && v14==v24){
    dir_swap=false; found=true; type=1;
  }
  if(v11==v22 && v12==v23 && v13==v24 && v14==v21){
    dir_swap=true; found=true; type=2;
  }
  if(v11==v22 && v12==v21 && v13==v24 && v14==v23){
    dir_swap=true; found=true; type=2;
  }
  if(v11==v23 && v12==v24 && v13==v21 && v14==v22){
    dir_swap=false; found=true; type=3;
  }
  if(v11==v24 && v12==v21 && v13==v22 && v14==v23){
    dir_swap=true; found=true; type=4;
  }
//  if(v11==v24 && v12==v23 && v13==v22 && v14==v21){
//    dir_swap=false; found=true; type=4;
  if(found==true){
    if(f1==1){ all_faces[i].left   =patch::to_string(j); e11=jmax[i]; e12=kmax[i];}
    if(f1==2){ all_faces[i].right  =patch::to_string(j); e11=jmax[i]; e12=kmax[i];}
    if(f1==3){ all_faces[i].bottom =patch::to_string(j); e11=imax[i]; e12=kmax[i];}
    if(f1==4){ all_faces[i].top    =patch::to_string(j); e11=imax[i]; e12=kmax[i];}
    if(f1==5){ all_faces[i].front  =patch::to_string(j); e11=imax[i]; e12=jmax[i];}
    if(f1==6){ all_faces[i].back   =patch::to_string(j); e11=imax[i]; e12=jmax[i];}
    if(f2==1){ all_faces[j].left   =patch::to_string(i); e21=jmax[j]; e22=kmax[j];}
    if(f2==2){ all_faces[j].right  =patch::to_string(i); e21=jmax[j]; e22=kmax[j];}
    if(f2==3){ all_faces[j].bottom =patch::to_string(i); e21=imax[j]; e22=kmax[j];}
    if(f2==4){ all_faces[j].top    =patch::to_string(i); e21=imax[j]; e22=kmax[j];}
    if(f2==5){ all_faces[j].front  =patch::to_string(i); e21=imax[j]; e22=jmax[j];}
    if(f2==6){ all_faces[j].back   =patch::to_string(i); e21=imax[j]; e22=jmax[j];}

    if(type==2){swap(s21,e21);type=1;}
    if(type==3){swap(s21,e21); swap(s22,e22);}
    if(type==4){swap(s22,e22);}
    if(mpi_class[i][f1]<0 && mpi_class[j][f2]<0){
      mpi_class[i][f1]=0;
      mpi_class[j][f2]=1;
    }
    write_to_mapfile(map,i,f1,s11,e11,s12,e12,j,f2,s21,e21,s22,e22,dir_swap);
  }
  
}

int main()
{
    std::ofstream block;
    block.open("blocking.dat");
    read_compute handler;
    handler.read_grid_names(block);
    block << "imx 	 jmx 	 kmx \n";
    for(int i =0; i< handler.tot_grids; i++)
    {
        cout << "Wait!! Reading grid "<< i+1 << endl;
        handler.compute_grid_medians(handler.grid_names[i], block);
    }
    cout << "computing adjacent blocks" << endl;
    handler.compute_interfaces();
    handler.write_to_file();
    block.close();
    //cout << handler.all_faces[0].left << handler.all_faces[0].right << endl;
    //cout << handler.all_vertices.size();
    //cout << handler.grid_names[1] << handler.tot_grids << endl;
    return 0;
}
