#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <fstream>

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
    std::vector <std::string> grid_names;
    std::vector <vertices> all_vertices;
    std::vector <faces> all_faces;
    void read_grid_names();//read file_names to get grids
    void compute_grid_medians(std::string str);
    void compute_interfaces();
    void write_to_file();
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
    output_file << i << "  " << grid_names[i]<<"  "<<"bc_"+patch::to_string(i)+".md" <<"  "<<all_faces[i].left << "  " << all_faces[i].right << "  "<< all_faces[i].bottom << "  ";
    output_file << all_faces[i].top << "  "<< all_faces[i].front << "  "<< all_faces[i].back << "  ";
    output_file << "\n";
    }
    output_file.close();

}

void read_compute::compute_interfaces()
{

    vertices buf1,buf2;
    coord c_buf1,c_buf2;
    for(int i =0; i< tot_grids ; i++)
    {
        /*
        left face   a -> c -> e -> g
        right face  b -> d -> f -> h
        back face   a -> c -> b -> d
        front face  e -> g -> f -> h
        bottom face a -> b -> e -> f
        top face c -> d -> g -> h
        */
        // left face check
        if(all_faces[i].left == "USER")
        {
           //cout << "enterd 1";
            buf1 = all_vertices[i];
            for(int j = 0; j< tot_grids; j++)
            {
                buf2 = all_vertices[j];
                c_buf1 = buf1.a;
                c_buf2 = buf2.b;
                //cout << c_buf1.x << " " << c_buf2.x << endl;
                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
                {

                    continue;
                }
                c_buf1 = buf1.c;
                c_buf2 = buf2.d;
                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
                {
                //cout << "enterd 10";
                    continue;
                }
                c_buf1 = buf1.e;
                c_buf2 = buf2.f;
                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
                {
                    continue;
                }
                c_buf1 = buf1.g;
                c_buf2 = buf2.h;
                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
                {
                    continue;
                }
                all_faces[i].left =patch::to_string(j);
                all_faces[j].right =patch::to_string(i);
                break;
            }
        }

        // top face check
        if(all_faces[i].bottom == "USER")
        {
            buf1 = all_vertices[i];
            for(int j = 0; j< tot_grids; j++)
            {
                buf2 = all_vertices[j];
                c_buf1 = buf1.a;
                c_buf2 = buf2.c;
                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
                {
                    continue;
                }
                c_buf1 = buf1.b;
                c_buf2 = buf2.d;
                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
                {
                    continue;
                }
                c_buf1 = buf1.e;
                c_buf2 = buf2.g;
                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
                {
                    continue;
                }
                c_buf1 = buf1.f;
                c_buf2 = buf2.h;
                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
                {
                    continue;
                }
                all_faces[i].bottom =patch::to_string(j);
                all_faces[j].top =patch::to_string(i);
                break;
            }
        }

        // front face check
        if(all_faces[i].front == "USER")
        {
            buf1 = all_vertices[i];
            for(int j = 0; j< tot_grids; j++)
            {
                buf2 = all_vertices[j];
                c_buf1 = buf1.e;
                c_buf2 = buf2.a;
                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
                {
                    continue;
                }
                c_buf1 = buf1.g;
                c_buf2 = buf2.c;
                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
                {
                    continue;
                }
                c_buf1 = buf1.f;
                c_buf2 = buf2.b;
                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
                {
                    continue;
                }
                c_buf1 = buf1.h;
                c_buf2 = buf2.d;
                if(c_buf1.x  != c_buf2.x || c_buf1.y  != c_buf2.y || c_buf1.z  != c_buf2.z)
                {
                    continue;
                }
                all_faces[i].front =patch::to_string(j);
                all_faces[j].back =patch::to_string(i);
                break;
            }
        }

    }

}
void read_compute::read_grid_names()
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
}

void read_compute::compute_grid_medians(std::string str)
{
    char *cr = &str[0u];
    std::ifstream grid;
    grid.open(cr);
    int imx =51,jmx =51, kmx =51;
    coord a,b,c,d,e,f,g,h,buf;
    medians buf_med;
    grid >> imx >> jmx >> kmx;
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
    grid.close();
}
std::string read_compute::layout_comment(string str)
{
return "## " + str + "\n";
}

int main()
{
    read_compute handler;
    handler.read_grid_names();
    for(int i =0; i< handler.tot_grids; i++)
    {
        cout << "Wait!! Reading grid "<< i+1 << endl;
        handler.compute_grid_medians(handler.grid_names[i]);
    }
    cout << "computing adjacent blocks" << endl;
    handler.compute_interfaces();
    handler.write_to_file();
    //cout << handler.all_faces[0].left << handler.all_faces[0].right << endl;
    //cout << handler.all_vertices.size();
    //cout << handler.grid_names[1] << handler.tot_grids << endl;
    return 0;
}
