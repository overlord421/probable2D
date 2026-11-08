#include <gapped2d_scattering.hh>
#include <algorithm>
#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <memory>
#include <numeric>
#include <sstream>
using namespace probable;
constexpr size_t NC=7,NH=60;
const std::array<double,NC> cuts{.02,.05,.1,.2,.5,1,2};
using Four=std::array<double,4>;
using Axes=std::array<Four,2>;
using Moments=std::array<Axes,NC>;
double logistic(double x){if(x>50)return std::exp(-x);if(x< -50)return 1-std::exp(x);return 1/(1+std::exp(x));}
double softplus(double x){return std::max(x,0.)+std::log1p(std::exp(-std::fabs(x)));}
struct Mode{std::string name;double Q,Do,nb;};
struct Model{
 Material m;double T,kT,rho,Da,sound,mu,nfl,density_check;std::vector<Mode> modes;
 std::vector<double> us,surfaces,ws;
 Model(double ncm,const std::string& modefile):m(1.285*consts::me,.125*consts::me,units::eV,.3,.3,100,300*units::K,300*units::K,'x',true,ncm*1e4/(units::m*units::m)),T(300*units::K),kT(consts::kB*T),rho(1.4e-3*units::g/(units::m*units::m)),Da(4.9*units::eV),sound(6800*units::m/units::s){
  std::ifstream input(modefile);if(!input)throw std::runtime_error("Cannot read mode file");
  std::string line;while(std::getline(input,line)){
   if(line.empty()||line[0]=='#')continue;
   std::istringstream row(line);std::string name;double energy,deformation;
   if(!(row>>name>>energy>>deformation)||!(energy>0)||!(deformation>0))throw std::runtime_error("Invalid mode row");
   double quantum=energy*1e-3*units::eV;
   modes.push_back({name,quantum,deformation*units::eV/units::m,1/std::expm1(quantum/kT)});
  }
  if(modes.empty())throw std::runtime_error("Empty optical mode set");
  const int N=8192;double hi=1.5*units::eV,h=hi/N;
  for(int j=0;j<=N;++j){double u=j*h;us.push_back(u);surfaces.push_back(surface(m.Delta+u));ws.push_back((j==0||j==N?1:j%2?4:2)*h/3);}
  double lo=m.Delta-60*kT,up=m.Delta+hi;
  for(int j=0;j<64;++j){double mid=(lo+up)/2;if(density(mid,false)<m.carrier_density_2d)lo=mid;else up=mid;}
  mu=(lo+up)/2;density_check=density(mu,false);nfl=density(mu,true);
 }
 double surface(double E)const{if(E<m.Delta)return 0;return 4*std::sqrt(m.mx*m.my*E/m.Delta)*math::comp_ellint_1(std::sqrt((E-m.Delta)/(2*E)));}
 double f(double E)const{return logistic((E-mu)/kT);}
 double density(double chem,bool fluct)const{
  double sum=0;for(size_t j=0;j<us.size();++j){double occ=logistic((m.Delta+us[j]-chem)/kT);sum+=ws[j]*surfaces[j]*occ*(fluct?1-occ:1);}
  return m.spin_degeneracy/std::pow(2*math::pi*consts::hbar,2)*sum;
 }
 std::vector<double> bare(double E)const{
  std::vector<double> result(1+2*modes.size());
  result[0]=2*Da*Da*kT/(math::pi*std::pow(consts::hbar,3)*rho*sound*sound)*surface(E)/4;
  for(size_t j=0;j<modes.size();++j){const auto& mode=modes[j];
   double op=mode.Do*mode.Do/(math::pi*std::pow(consts::hbar,2)*rho*(mode.Q/consts::hbar));
   result[1+2*j]=E-mode.Q>=m.Delta?op*surface(E-mode.Q)/4*(mode.nb+1):0;
   result[2+2*j]=op*surface(E+mode.Q)/4*mode.nb;
  }
  return result;
 }
 std::vector<double> rates(double E)const{
  auto b=bare(E);
  for(size_t j=0;j<modes.size();++j)for(size_t direction=0;direction<2;++direction){
   size_t channel=1+2*j+direction;double shift=direction==0?modes[j].Q:-modes[j].Q;
   if(b[channel]>0)b[channel]*=std::exp(softplus((mu-E)/kT)-softplus((mu-(E-shift))/kT));
  }
  return b;
 }
 Vec2 sample(const Scattering& angular)const{
  double emax=m.Delta+std::max(0.,mu-m.Delta)+32*kT;
  double flo=f(emax),fhi=f(m.Delta),Smax=surface(emax);
  while(true){double occ=flo+(fhi-flo)*uniform();double E=mu+kT*std::log((1-occ)/occ);
   if(uniform()*Smax<=surface(E)){
    Vec2 p{std::sqrt(m.mx*(E-m.Delta)*(E+m.Delta)/m.Delta),0};return angular.scatter(p);
   }
  }
 }
};
void add(Axes& out,Vec2 v,Vec2 v0,double q,double q0,double dt){
 std::array<double,2>a{v.x,v.y},b{v0.x,v0.y};
 for(size_t ax=0;ax<2;++ax){double z=a[ax]*b[ax]*dt;out[ax][0]+=z;out[ax][1]+=q*z;out[ax][2]+=q0*z;out[ax][3]+=q*q0*z;}
}
void print_moments(const Moments&mom,size_t N){
 std::cout<<"[";for(size_t c=0;c<NC;++c){if(c)std::cout<<",";std::cout<<"[";for(size_t ax=0;ax<2;++ax){if(ax)std::cout<<",";std::cout<<"[";for(size_t j=0;j<4;++j){if(j)std::cout<<",";std::cout<<mom[c][ax][j]/N;}std::cout<<"]";}std::cout<<"]";}std::cout<<"]";
}
int main(int argc,char**argv){
 if(argc!=4)return 2;double ncm=std::stod(argv[1]);size_t N=std::stoul(argv[2]);Model model(ncm,argv[3]);
 AcousticScattering ac(model.m,model.rho,model.Da,model.sound);
 std::vector<std::unique_ptr<Scattering>> optical;
 std::vector<Scattering*> channels{&ac};
 for(const auto& mode:model.modes){
  optical.emplace_back(new OpticalEmissionScattering(model.m,model.rho,mode.Do,mode.Q));channels.push_back(optical.back().get());
  optical.emplace_back(new OpticalAbsorptionScattering(model.m,model.rho,mode.Do,mode.Q));channels.push_back(optical.back().get());
 }
 Moments full{},rb{};std::vector<unsigned long long> events(channels.size());
 std::array<double,NH> residence{},initial_hist{},final_hist{};double max_energy_error=0,max_db=0,rate_change=0;
 double histmax=.6*units::eV;
 auto bin=[&](double E){return std::min(NH-1,size_t(std::max(0.,E-model.m.Delta)/histmax*NH));};
 auto start=std::chrono::steady_clock::now();
 for(double u:{0.,.005,.0417864,.0417867,.0544644,.0544647,.12,.2,.4}){
  double E=model.m.Delta+u*units::eV;auto a=model.bare(E);
  for(size_t j=0;j<model.modes.size();++j){double Q=model.modes[j].Q;auto b=model.bare(E+Q);
   double lhs=model.surface(E)*model.f(E)*a[2+2*j]*(1-model.f(E+Q));
   double rhs=model.surface(E+Q)*model.f(E+Q)*b[1+2*j]*(1-model.f(E));
   max_db=std::max(max_db,std::fabs(lhs-rhs)/std::max(lhs+rhs,1e-300));
  }
  Vec2 p{std::sqrt(model.m.mx*(E-model.m.Delta)*(E+model.m.Delta)/model.m.Delta),0};
  if(u<=.2)rate_change=std::max(rate_change,std::fabs(ac.rate(p,0)/a[0]-1));
 }
 for(size_t j=0;j<N;++j){
  Vec2 p=model.sample(ac),v0=model.m.velocity(p);double E0=model.m.energy(p),q0=E0-model.mu;
  initial_hist[bin(E0)]++;
  auto rates0=model.rates(E0);double gamma0=std::accumulate(rates0.begin(),rates0.end(),0.);
  for(size_t c=0;c<NC;++c){double duration=-std::expm1(-gamma0*cuts[c]*units::ps)/gamma0;add(rb[c],v0,v0,q0,q0,duration);}
  Moments path{};double time=0;size_t c=0;
  while(c<NC){
   double E=model.m.energy(p),q=E-model.mu;Vec2 v=model.m.velocity(p);auto rates=model.rates(E);
   double gamma=std::accumulate(rates.begin(),rates.end(),0.);if(!(gamma>0)||!std::isfinite(gamma))throw std::runtime_error("Bad response rate");
   double event=time-std::log(uniform())/gamma;
   while(c<NC&&cuts[c]*units::ps<=event){double stop=cuts[c]*units::ps;Axes segment{};add(segment,v,v0,q,q0,stop-time);
    for(size_t a=0;a<2;++a)for(size_t b=0;b<4;++b)path[c][a][b]+=segment[a][b];
    residence[bin(E)]+=stop-time;time=stop;
    if(c+1<NC)path[c+1]=path[c];++c;
   }
   if(c==NC)break;
   Axes segment{};add(segment,v,v0,q,q0,event-time);for(size_t a=0;a<2;++a)for(size_t b=0;b<4;++b)path[c][a][b]+=segment[a][b];
   residence[bin(E)]+=event-time;time=event;
   double select=uniform()*gamma,cumulative=rates[0];size_t ch=0;while(ch+1<rates.size()&&select>=cumulative)cumulative+=rates[++ch];
   Vec2 next=channels[ch]->scatter(p);double target=E-channels[ch]->energy;
   max_energy_error=std::max(max_energy_error,std::fabs(model.m.energy(next)-target)/units::eV);p=next;events[ch]++;
  }
  for(size_t ci=0;ci<NC;++ci)for(size_t ax=0;ax<2;++ax)for(size_t k=0;k<4;++k)full[ci][ax][k]+=path[ci][ax][k];
  final_hist[bin(model.m.energy(p))]++;
 }
 std::cout<<std::setprecision(17)<<"{\"n_cm2\":"<<ncm<<",\"N\":"<<N<<",\"T_K\":300,\"Delta_eV\":1"
 <<",\"mu_eV\":"<<model.mu/units::eV<<",\"density_ratio\":"<<model.density_check/model.m.carrier_density_2d
 <<",\"nfl_m2\":"<<model.nfl*units::m*units::m<<",\"density_m2\":"<<model.m.carrier_density_2d*units::m*units::m
 <<",\"sigma_prefactor_S\":"<<model.nfl*consts::e*consts::e/(consts::kB*model.T)/(units::C/units::V/units::s)
 <<",\"kappa_prefactor_W_mK\":"<<model.nfl/(consts::kB*model.T*model.T)/(units::J/units::s/units::K)/(.53e-9)
 <<",\"seebeck_prefactor_V_K\":"<<-1/(consts::e*model.T)/(units::V/units::K)
 <<",\"cutoffs_ps\":[";for(size_t c=0;c<NC;++c){if(c)std::cout<<",";std::cout<<cuts[c];}
 std::cout<<"],\"full_moments\":";print_moments(full,N);std::cout<<",\"conditional_moments\":";print_moments(rb,N);
 std::cout<<",\"events\":[";for(size_t j=0;j<events.size();++j){if(j)std::cout<<",";std::cout<<events[j];}
 std::cout<<"],\"max_energy_error_eV\":"<<max_energy_error
 <<",\"max_pair_balance_relative_residual\":"<<max_db<<",\"max_bare_rate_change_below_0p2eV\":"<<rate_change
 <<",\"histogram_max_eV\":0.6,\"residence_ps\":[";for(size_t j=0;j<NH;++j){if(j)std::cout<<",";std::cout<<residence[j];}
 std::cout<<"],\"initial_histogram\":[";for(size_t j=0;j<NH;++j){if(j)std::cout<<",";std::cout<<initial_hist[j];}
 std::cout<<"],\"final_histogram\":[";for(size_t j=0;j<NH;++j){if(j)std::cout<<",";std::cout<<final_hist[j];}
 std::cout<<"],\"elapsed_s\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
}
