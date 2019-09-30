#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <utility>

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

struct all_bc
{
    int imin;
    int imax;
    int jmin;
    int jmax;
    int kmin;
    int kmax;

};

class read_write
{
public:
    std::map <int,std::string> bc_values;
    std::map <int,std::string> bc_list;
    void read_bc_values();
    void read_bc_list();
    void read_layout();
    void write_bc_file(std::string, all_bc);
    int check_for_periodic_bc(std::ostream& periodic, int bc,  int count, int face );
    std::string parse_bc(std::string);

};
std::string read_write::parse_bc(std::string bc)
{
    std::string buf ="",net_string="";
    int bc_num;
    for(int i =0; i< bc.length(); i++)
    {
        if(bc[i] == '-')
        {

            bc_num =atoi(buf.c_str());
            net_string =net_string+ "- " + bc_values[bc_num] +"\n";
            buf = "";
        }
        else
            buf = buf+bc[i];
    }
    bc_num =atoi(buf.c_str());
    net_string =net_string+ "- " + bc_values[bc_num] +"\n\n";
    buf = "";
    return net_string;
}

void read_write::write_bc_file(std::string bc_file, all_bc buf)
{
    string str = "../bc/"+bc_file;
    char *c = &str[0u];
    std::fstream output_file (c,fstream::out);
    if(!output_file.is_open())
    {
        cout << "Error: Couldn't open output file" << endl;
        return;
    }
    output_file << "BOUNDARY CONDITIONS CONFIGURATION\n=================================\n\n";

    std::map<int,std::string>::iterator it;
    std::string bc;
    int bc_num;

    // imin
    bc_num = buf.imin;
    it = bc_list.find(bc_num);
    if(it == bc_list.end() && bc_num < 0)
    {
        cout << "ERROR: BC '"<<bc_num<< "' IN LIST NOT FOUND";
    }
    if(bc_num >= 0)
        bc_num = 1;
    output_file <<"# imn\n";
    bc = bc_list[bc_num];
    output_file << parse_bc(bc);

     // imax
    bc_num = buf.imax;
    it = bc_list.find(bc_num);
    if(it == bc_list.end() && bc_num < 0)
    {
        cout << "ERROR: BC '"<<bc_num<< "' IN LIST NOT FOUND";
    }
    if(bc_num >= 0)
        bc_num = 1;
    output_file <<"# imx\n";
    bc = bc_list[bc_num];
    output_file << parse_bc(bc);

         // jmin
    bc_num = buf.jmin;
    it = bc_list.find(bc_num);
    if(it == bc_list.end() && bc_num < 0)
    {
        cout << "ERROR: BC '"<<bc_num<< "' IN LIST NOT FOUND";
    }
    if(bc_num >= 0)
        bc_num = 1;
    output_file <<"# jmn\n";
    bc = bc_list[bc_num];
    output_file << parse_bc(bc);

         // jmax
    bc_num = buf.jmax;
    it = bc_list.find(bc_num);
    if(it == bc_list.end() && bc_num < 0)
    {
        cout << "ERROR: BC '"<<bc_num<< "' IN LIST NOT FOUND";
    }
    if(bc_num >= 0)
        bc_num = 1;
    output_file <<"# jmx\n";
    bc = bc_list[bc_num];
    output_file << parse_bc(bc);

             // kmin
    bc_num = buf.kmin;
    it = bc_list.find(bc_num);
    if(it == bc_list.end() && bc_num < 0)
    {
        cout << "ERROR: BC '"<<bc_num<< "' IN LIST NOT FOUND";
    }
    if(bc_num >= 0)
        bc_num = 1;
    output_file <<"# kmn\n";
    bc = bc_list[bc_num];
    output_file << parse_bc(bc);

         // kmax
    bc_num = buf.kmax;
    it = bc_list.find(bc_num);
    if(it == bc_list.end() && bc_num < 0)
    {
        cout << "ERROR: BC '"<<bc_num<< "' IN LIST NOT FOUND";
    }
    if(bc_num >= 0)
        bc_num = 1;
    output_file <<"# kmx\n";
    bc = bc_list[bc_num];
    output_file << parse_bc(bc);
    output_file << "FIN\n";



}

void read_write::read_layout()
{
    std::ofstream periodic;
    periodic.open("periodic.txt");
    periodic << "   B1   B2   F1   F2 "<<endl;
    fstream fin("layout.md");
    std::string line;
    int num,count=0,tot_proc,tot_ent;
    int proc_id;
    std::string imin,imax,jmin,jmax,kmin,kmax;
    //std::string imin="USER",imax="USER",jmin="USER",jmax="USER",kmin="USER",kmax="USER";
    std::string grid,bc_file;
    all_bc buf;
    while(getline(fin, line))
    {
        //the following line trims white space from the beginning of the string
        line.erase(line.begin(), find_if(line.begin(), line.end(), not1(ptr_fun<int, int>(isspace))));

        if(line[0] == '#') continue;//ignore lines starting with #
        count++;

        if(count > 2)
        {
            // boundary condition line
            stringstream(line) >> proc_id >> grid >> bc_file >> imin >> imax >> jmin >> jmax >> kmin >> kmax;
             if(imin == "USER" || imax == "USER" || jmin =="USER" ||jmax =="USER" || kmin == "USER" || kmax =="USER")
            {
               cout << "ERROR: ONE OR MORE USER BCs ARE YET TO BE GIVEN FOR PROCESS " << proc_id << endl;
               cout << "EROOR: NOT GENERATING BC FILE " << bc_file << endl;
               continue;
            }

            buf.imin = atoi(imin.c_str());
            buf.imax = atoi(imax.c_str());
            buf.jmin = atoi(jmin.c_str());
            buf.jmax = atoi(jmax.c_str());
            buf.kmin = atoi(kmin.c_str());
            buf.kmax = atoi(kmax.c_str());
            buf.imin = check_for_periodic_bc(periodic, buf.imin, count-3, 1);
            buf.jmin = check_for_periodic_bc(periodic, buf.jmin, count-3, 3);
            buf.kmin = check_for_periodic_bc(periodic, buf.kmin, count-3, 5);
            buf.imax = check_for_periodic_bc(periodic, buf.imax, count-3, 2);
            buf.jmax = check_for_periodic_bc(periodic, buf.jmax, count-3, 4);
            buf.kmax = check_for_periodic_bc(periodic, buf.kmax, count-3, 6);
            write_bc_file(bc_file,buf);

        }
        else if(count == 1)
            stringstream(line) >> tot_proc;
        else if(count == 2)
            stringstream(line) >> tot_ent;
    }
    //cout << tot_proc << " " << tot_ent <<endl;

}

int read_write::check_for_periodic_bc(std::ostream& periodic, int bc, int count, int face){
  // get the first digit
  std::string s = std::to_string(bc);
  int periodicF2;
  if(bc<0 && s.length()>2 && s.substr(1,1)==std::to_string(2)) {
    //write periodic bc control file
    if(face == 1) {periodicF2 = 2;}
    if(face == 2) {periodicF2 = 1;}
    if(face == 3) {periodicF2 = 4;}
    if(face == 4) {periodicF2 = 3;}
    if(face == 5) {periodicF2 = 6;}
    if(face == 6) {periodicF2 = 5;}
    periodic << " " << setw(4) << count << " " << setw(4) << s.substr(2) 
             << " " << setw(4) << face << " " << setw(4) << periodicF2 << endl;
    return -10;
    }
  else{ 
    return bc;
    }
  }

void read_write::read_bc_values()
{
    fstream fin("bc_values");
    std::string line,bc;
    int num;
    while(getline(fin, line))
    {
        //the following line trims white space from the beginning of the string
        line.erase(line.begin(), find_if(line.begin(), line.end(), not1(ptr_fun<int, int>(isspace))));

        if(line[0] == '#') continue;//ignore lines starting with #


        stringstream(line) >> num >> bc;
        bc_values[num] = bc;
        //cout << "Data: " << num << bc << endl;
    }
}

void read_write::read_bc_list()
{
    fstream fin("bc_list");
    std::string line,bc;
    int num;
    while(getline(fin, line))
    {
        //the following line trims white space from the beginning of the string
        line.erase(line.begin(), find_if(line.begin(), line.end(), not1(ptr_fun<int, int>(isspace))));

        if(line[0] == '#') continue;//ignore lines starting with #


        stringstream(line) >> num >> bc;
        bc_list[num] = bc;
        //cout << "Data: " << num << bc << endl;
    }
}


int main()
{
    read_write handler;
    handler.read_bc_values();
    handler.read_bc_list();
    handler.read_layout();
    return 0;
}
