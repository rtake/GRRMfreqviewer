//  Written by 計算化学.com 管理人 2/1/2017
//  This program converts GRRM's Freq log file to Gaussian's format.
//
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdio>
using namespace std;
int genso(string atom);    //prototype declaration
int main(int argc, char *argv[]){
//  Input Error
if(argc == 1){
cout << "\n *************************\n ERROR -- no argumanet --\n ./a.out [file_name]\n *************************\n\n Copyright (c) 2017\n 計算化学.com  All Rights Reserved.\n" << endl;
exit(1);
}
//  Reading file name
ifstream ifs(argv[1]);
string line;
int i(0), k(0), Geometry_startline(0), Initial_Structure(0), Gradient_Startline(0), Gradient_Endline(0);
long Freq_startline[9999];
while (getline(ifs,line)){
i++;
std::string::size_type Geometry_Find = line.find("Geometry");
std::string::size_type Freq_Find = line.find("Freq.");
std::string::size_type Initial_Find = line.find("INPUT ORIENTATION");
std::string::size_type Gradient_Find = line.find("GRADIENT VECTOR");
std::string::size_type Hessian_Find = line.find("HESSIAN MATRIX");
if (Initial_Find != std::string::npos ){
Initial_Structure = i;
}if (Gradient_Find != std::string::npos ){
Gradient_Startline = i;
}if (Hessian_Find != std::string::npos ){
Gradient_Endline = i;
}if (Geometry_Find != std::string::npos ){
Geometry_startline = i;
}if (Freq_Find != std::string::npos ) {
k++;
Freq_startline[k] = i;
}
}
int atom_number = Freq_startline[1] - Geometry_startline -3;
//writing outputfile
char filename_without_com[50];
int file_length (0);
for(int i (0); argv[1][i] !='\0'; ++i) file_length +=1;
for(int i = 0; i <= file_length-4; i++) filename_without_com[i] = argv[1][i];
filename_without_com[file_length-4] = '\0';
stringstream outputfile_name;
outputfile_name << filename_without_com << "_convert.log";
ofstream ofs(outputfile_name.str());
ofs <<  "----------------------------------------------------------------------\n#P\n----------------------------------------------------------------------\n\n" << "  and normal coordinates:" << endl;
ifs.close();
ifs.open(argv[1]);
i=0, k=1;
int l(0);
double vibration[99999][3]={};
string element[99999];
while (getline(ifs,line)){
i++;
if(i == Freq_startline[k]-1){
int AA(0), BB(0), CC(0);
sscanf(line.data(),"%d %d %d", &AA, &BB, &CC);
ofs << "                      " << AA +1 << "                      " << BB +1 << "                      " << CC +1 <<endl;
ofs << "                      A                      A                      A" << endl;
}if(i == Freq_startline[k]){
double Frequencies[3]={};
sscanf(line.substr(10,47).data(),"%lf %lf %lf", &Frequencies[0], &Frequencies[1], &Frequencies[2]);
ofs << "Frequencies --     " << Frequencies[0] << "                " << Frequencies[1] << "               " << Frequencies[2] << endl;
ofs << "  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z"<<endl;
}if(i > Freq_startline[k] + 1 && i < (Freq_startline[k] + 1 + atom_number * 3) ){
l++;
double a(0), b(0), c(0);
sscanf(line.substr(10,47).data(),"%lf %lf %lf", &a, &b, &c);
vibration[l][0] = a;
vibration[l][1] = b;
vibration[l][2] = c;
element[l] = line.substr(0,2);
}if(i == (Freq_startline[k] + 1 + atom_number * 3)){
for(int j=0; j<=atom_number-1; j++){
//writing element symbol
if(j+1<10) ofs<< "     " << j+1;
else ofs<< "    " << j+1;
if (genso(element[j*3+1])<10) ofs << "   " << genso(element[j*3+1]) << "    ";
else ofs << "  " << genso(element[j*3+1]) << "    ";
//1つ目
if(vibration[j*3+1][0]>=0) ofs << " ";
ofs << fixed << setprecision(2) << vibration[j*3+1][0] << "  ";
if(vibration[j*3+2][0]>=0) ofs << " ";
ofs << fixed << setprecision(2) << vibration[j*3+2][0] << "  ";
if(vibration[j*3+3][0]>0) ofs << " ";
ofs << fixed << setprecision(2) << vibration[j*3+3][0] << "\t";
//2つ目
if(vibration[j*3+1][1]>=0) ofs << " ";
ofs << fixed << setprecision(2) << vibration[j*3+1][1] << "  ";
if(vibration[j*3+2][1]>=0) ofs << " ";
ofs << fixed << setprecision(2) << vibration[j*3+2][1] << "  ";
if(vibration[j*3+3][1]>0) ofs << " ";
ofs << fixed << setprecision(2) << vibration[j*3+3][1] << "\t";
//3つ目
if(vibration[j*3+1][2]>=0) ofs << " ";
ofs << fixed << setprecision(2) << vibration[j*3+1][2] << "  ";
if(vibration[j*3+2][2]>=0) ofs << " ";
ofs << fixed << setprecision(2) << vibration[j*3+2][2] << "  ";
if(vibration[j*3+3][2]>0) ofs << " ";
ofs << fixed << setprecision(2) << vibration[j*3+3][2] << endl;
}
k++;
}
}
ofs << "\n***** Axes restored to original set *****\n" << " -------------------------------------------------------------------\n" << " Center     Atomic                   Forces (Hartrees/Bohr)\n" << " Number     Number              X              Y              Z\n" << " -------------------------------------------------------------------\n";
// writing gradient
ifs.close();
ifs.open(argv[1]);
double gradient[99999];
i=0, k=0;
while (getline(ifs,line)){
i++;
if(i > Gradient_Startline && i < Gradient_Endline){
k++;
double a(0);
sscanf(line.data(),"%lf", &a);
gradient[k] = a;
}
}
for(int j=0; j<=atom_number-1; j++){
//writing element symbol
if(j+1<10) ofs<< "     " << j+1;
else ofs<< "    " << j+1;
//atomic number
if (genso(element[j*3+1])<10) ofs << "        " << genso(element[j*3+1]) << "         ";
else ofs << "       " << genso(element[j*3+1]) << "         ";
//1つ目
if(gradient[j*3+1]>=0) ofs << " ";
ofs << fixed << setprecision(10) << gradient[j*3+1] << "   ";
if(gradient[j*3+2]>=0) ofs << " ";
ofs << fixed << setprecision(10) << gradient[j*3+2] << "   ";
if(gradient[j*3+3]>0) ofs << " ";
ofs << fixed << setprecision(10) << gradient[j*3+3] << endl;
}
//initial coordinate
ifs.close();
ifs.open(argv[1]);
i=0;
stringstream Coordinates;
while (getline(ifs,line)){
i++;
if(i == Initial_Structure){
ofs<<" Test job not archived.\n"<<" 1\\1\\G\\OPT\\R\\C\\S\n"<<" 5\\0\\\\#P\\"<<endl;
Coordinates<< " q\\\\t\\\\1,1";
}if(i == Initial_Structure+1){
double a(0), b(0), c(0);
sscanf(line.substr(3,54).data(),"%lf\t%lf\t%lf",&a,&b,&c);
if(line.substr(1,1) == " ") Coordinates<<"\\"<<line.substr(0,1)<<","<<a<<","<<b<<","<<c;
else Coordinates<<"\\"<<line.substr(0,2)<<","<<a<<","<<b<<","<<c;
}if(i > Initial_Structure+1 && i < Initial_Structure + atom_number){
double a(0), b(0), c(0);
sscanf(line.substr(3,54).data(),"%lf\t%lf\t%lf",&a,&b,&c);
if(line.substr(1,1) == " ") Coordinates<<"\\"<<line.substr(0,1)<<","<<a<<","<<b<<","<<c;
else Coordinates<<"\\"<<line.substr(0,2)<<","<<a<<","<<b<<","<<c;
}if(i == Initial_Structure + atom_number){
double a(0), b(0), c(0);
sscanf(line.substr(3,54).data(),"%lf\t%lf\t%lf",&a,&b,&c);
if(line.substr(1,1) == " ") Coordinates<<"\\"<<line.substr(0,1)<<","<<a<<","<<b<<","<<c;
else Coordinates<<"\\"<<line.substr(0,2)<<","<<a<<","<<b<<","<<c;
}
}
Coordinates<<"\\\\Version=Fujitsu-XTC-G16RevB.01\\";
for(int j=0; j<=(Coordinates.str().length()/70) ; j++){
ofs << " " << Coordinates.str().substr(70*j,70)<<endl;
}
ofs<<"[X(C20H37O2)]\\\\\\@"<<endl<<endl;
ofs<<"Normal termination of Gaussian 16"<<endl;
stringstream buff;
buff << "open /Applications/gv/gview.app " << outputfile_name.str() << " &\n";
//cout << buff.str().c_str() << endl;
system(buff.str().c_str());
ofs.close();
}
//End of Main Functional
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function converting element name to atomic nuimber
int genso (string atom){
string elements[128]={"Error","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Al","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Ti","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Md","Bo","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cp","Uut","Uuq","Uup","Uuh","Uus","Uus","Uuo"};
int j = 1;
if(atom.substr(1,1)==" "){
for(int i(1); i<=118; i++){
if(strcmp(elements[i].c_str(),(atom.substr(0,1)).c_str())==0){
j = i;break;
}
else continue;}
}else{
for(int i(1); i<=118; i++){
if(strcmp(elements[i].c_str(),(atom.substr(0,2)).c_str())==0){
j = i;break;
}
else continue;}}
return j;
}
