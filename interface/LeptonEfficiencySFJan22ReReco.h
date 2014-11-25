#ifndef LeptonEfficiencySF_h
#define LeptonEfficiencySF_h

// cf.
// https://twiki.cern.ch/twiki/bin/view/Main/EGammaScaleFactors2012#2012_8_TeV_Jan22_Re_recoed_data
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs#22Jan2013_ReReco_of_2012_data_re
//
class LeptonEfficiencySF
{
 public:

  //
  LeptonEfficiencySF() { }

  //
  ~LeptonEfficiencySF() {}

  //
  std::pair<float,float> getLeptonEfficiency(float pt, float eta, int id, std::string wp)
    {
      eta=fabs(eta);
      id=abs(id);
      
      std::pair<float,float> eff(1.0,0.04);
      switch(id){
      case 11:
	{
        if(wp=="loose")
	    {
            if(fabs(eta)<0.8)
            {
                if(pt<30)      { eff.first=0.988;  eff.second=0.002; }
                else if(pt<40) { eff.first=1.002;  eff.second=0.001; }
                else if(pt<50) { eff.first=1.005;  eff.second=0.001; }
                else           { eff.first=1.005;  eff.second=0.001; }

            }
            else if(fabs(eta)<1.442)
            {
                if(pt<30)      { eff.first=0.965;  eff.second=0.003; }
                else if(pt<40) { eff.first=0.985;  eff.second=0.001; }
                else if(pt<50) { eff.first=0.989;  eff.second=0.001; }
                else           { eff.first=0.989;  eff.second=0.002; }
            }
            else if(fabs(eta)<1.556)
            {
                if(pt<30)      { eff.first=0.990;  eff.second=0.011; }
                else if(pt<40) { eff.first=0.966;  eff.second=0.005; }
                else if(pt<50) { eff.first=0.971;  eff.second=0.004; }
                else           { eff.first=0.980;  eff.second=0.008; }
            }
            else if(fabs(eta)<2.0)
            {
                if(pt<30)      { eff.first=0.953;  eff.second=0.005; }
                else if(pt<40) { eff.first=0.980;  eff.second=0.003; }
                else if(pt<50) { eff.first=0.999;  eff.second=0.002; }
                else           { eff.first=1.004;  eff.second=0.0025; }
            }
            else
            {
                if(pt<30)      { eff.first=1.017;  eff.second=0.005; }
                else if(pt<40) { eff.first=1.019;  eff.second=0.003; }
                else if(pt<50) { eff.first=1.019;  eff.second=0.002; }
                else           { eff.first=1.023;  eff.second=0.004; }
            }
        }
        
        if(wp=="medium")
        {
            if(fabs(eta)<0.8)
            {
                if(pt<30)      { eff.first=0.986;  eff.second=0.0015; }
                else if(pt<40) { eff.first=1.002;  eff.second=0.001; }
                else if(pt<50) { eff.first=1.005;  eff.second=0.001; }
                else           { eff.first=1.004;  eff.second=0.001; }
            }
            else if(fabs(eta)<1.442)
            {
                if(pt<30)      { eff.first=0.959;  eff.second=0.003; }
                else if(pt<40) { eff.first=0.980;  eff.second=0.001; }
                else if(pt<50) { eff.first=0.988;  eff.second=0.001; }
                else           { eff.first=0.988;  eff.second=0.002; }
            }
            else if(fabs(eta)<1.556)
            {
                if(pt<30)      { eff.first=0.967;  eff.second=0.010; }
                else if(pt<40) { eff.first=0.950;  eff.second=0.0065; }
                else if(pt<50) { eff.first=0.958;  eff.second=0.005; }
                else           { eff.first=0.966;  eff.second=0.009; }
            }
            else if(fabs(eta)<2.0)
            {
                if(pt<30)      { eff.first=0.941;  eff.second=0.005; }
                else if(pt<40) { eff.first=0.967;  eff.second=0.003; }
                else if(pt<50) { eff.first=0.992;  eff.second=0.002; }
                else           { eff.first=1.000;  eff.second=0.003; }
            }
            else
            {
                if(pt<30)      { eff.first=1.020;  eff.second=0.003; }
                else if(pt<40) { eff.first=1.021;  eff.second=0.003; }
                else if(pt<50) { eff.first=1.019;  eff.second=0.002; }
                else           { eff.first=1.022;  eff.second=0.004; }
            }
        }
	}
	break;
      case 13:
	{
	  if(wp=="loose")
	    {

            if(fabs(eta)<0.9)
            {
                if(pt<20)	     { eff.first=0.975506202184*0.963643517924;  eff.second=0.00553198576496; }
                else if(pt<25)       { eff.first=0.9862024644*0.988643237368;   eff.second=0.00195537767007; }
                else if(pt<30)       { eff.first=0.993370286278*0.999383313732;   eff.second=0.000906009608215; }
                else if(pt<35)       { eff.first=0.996768819346*0.998655709693;   eff.second=0.000549384532056; }
                else if(pt<40)       { eff.first=0.997442382829*0.998344487569;   eff.second=0.000346171326837; }
                else if(pt<50)       { eff.first=0.998198010419*0.998169667151;   eff.second=0.000167634519736; }
                else if(pt<60)       { eff.first=0.995041805799*0.99906381197;   eff.second=0.00039408352652; }
                else if(pt<90)       { eff.first=0.987628409023*1.00043605512;   eff.second=0.000710617117054; }
                else if(pt<140)      { eff.first=0.998926292328*1.00080008919;   eff.second=0.00242005780083; }
                else if(pt<300)      { eff.first=1.00230988809*1.00188008471;   eff.second=0.0126337604302; }
            }
            else if(fabs(eta)<1.2)
            {
                if(pt<20)	     { eff.first=1.00452783959*0.963850533573;  eff.second=0.00462179564672; }
                else if(pt<25)       { eff.first=0.999234569533*0.987899669049;   eff.second=0.00280811707302; }
                else if(pt<30)       { eff.first=0.998902784637*1.00147584363;   eff.second=0.00154831047279; }
                else if(pt<35)       { eff.first=0.997492892443*1.00213211635;   eff.second=0.00103035387795; }
                else if(pt<40)       { eff.first=0.998633235305*1.00105023257;   eff.second=0.000612766325709; }
                else if(pt<50)       { eff.first=0.99810724692*1.00004986605;   eff.second=0.00028329810884; }
                else if(pt<60)       { eff.first=0.997018342883*0.999964173777;   eff.second=0.000691478840647; }
                else if(pt<90)       { eff.first=0.989447339964*1.00101227493;   eff.second=0.001313332479; }
                else if(pt<140)      { eff.first=1.00105122689*1.00153364615;   eff.second=0.00350893654878; }
                else if(pt<300)      { eff.first=1.00430604383*0.998248497453;   eff.second=0.0108062348387; }
            }
            else if(fabs(eta)<2.1)
            {
                if(pt<20)	     { eff.first=1.00220514615*0.977631976574;  eff.second=0.00178024447192; }
                else if(pt<25)       { eff.first=1.0031908843*0.994508796984;   eff.second=0.0011891203869; }
                else if(pt<30)       { eff.first=0.998140970618*1.00207203892;   eff.second=0.000796368030622; }
                else if(pt<35)       { eff.first=0.998022947197*1.00289300399;   eff.second=0.000583381354257; }
                else if(pt<40)       { eff.first=0.998220167836*1.00184531204;   eff.second=0.000391211404999; }
                else if(pt<50)       { eff.first=0.998350291884*1.00000294414;   eff.second=0.000186035995102; }
                else if(pt<60)       { eff.first=0.994845625243*1.00011257623;   eff.second=0.00048651408496; }
                else if(pt<90)       { eff.first=0.984052544696*1.00037169484;   eff.second=0.00101161306492; }
                else if(pt<140)      { eff.first=1.00429862205*0.999769938901;   eff.second=0.0019985241959; }
                else if(pt<300)      { eff.first=0.98537402128*0.997082888406;   eff.second=0.0200267470732; }
            }
            else if(fabs(eta)<2.4)
            {
                if(pt<20)            { eff.first=1.00762887238*1.06688690031;  eff.second=0.00432394972324; }
                else if(pt<25)       { eff.first=1.00464891352*1.05397476812;   eff.second=0.00339452951027; }
                else if(pt<30)       { eff.first=1.00156206018*1.04187105417;   eff.second=0.00201582420415; }
                else if(pt<35)       { eff.first=1.00220667892*1.02812749722;   eff.second=0.00146850313709; }
                else if(pt<40)       { eff.first=0.998205637346*1.02036051109;   eff.second=0.00109927230931; }
                else if(pt<50)       { eff.first=0.998660592368*1.00934163025;   eff.second=0.000633160288086; }
                else if(pt<60)       { eff.first=0.992784180151*1.00647775397;   eff.second=0.00164121434978; }
                else if(pt<90)       { eff.first=0.972010302563*1.00471548378;   eff.second=0.00361440532758; }
                else if(pt<140)      { eff.first=1.01574040111*0.999605862868;   eff.second=0.00997644433652; }
                else if(pt<300)      { eff.first=0.948823376501*1.01096438407;   eff.second=0.100546441794; }
            }

	    }
	  if(wp=="tight"){
          
          if(fabs(eta)<0.9)
            {
                if(pt<20)            { eff.first=0.970274891725*0.958554337961;  eff.second=0.00567231311004; }
                else if(pt<25)       { eff.first=0.988864710503*0.988192642644;   eff.second=0.00219145424653; }
                else if(pt<30)       { eff.first=0.992338228944*0.999494808945;   eff.second=0.00109858582064; }
                else if(pt<35)       { eff.first=0.993283243571*0.998766213813;   eff.second=0.000725620185173; }
                else if(pt<40)       { eff.first=0.993661904524*0.998527321041;   eff.second=0.00051603639698; }
                else if(pt<50)       { eff.first=0.99239188725*0.998498730357;   eff.second=0.000300488942088; }
                else if(pt<60)       { eff.first=0.991189978664*0.999247195464;   eff.second=0.000703641521117; }
                else if(pt<90)       { eff.first=0.989416795661*1.00046571514;   eff.second=0.00108987827443; }
                else if(pt<140)      { eff.first=1.00374898772*1.00074782308;   eff.second=0.00324108730088; }
                else if(pt<300)      { eff.first=1.01850255404*1.00090397445;   eff.second=0.0174696418775; }
            }
          else if(fabs(eta)<1.2)
          {
              	if(pt<20)	     { eff.first=1.00173130936*0.96628645395;  eff.second=0.00742746410077; }
              	else if(pt<25)       { eff.first=0.993946645189*0.989546384934;   eff.second=0.00319897078553; }
              	else if(pt<30)       { eff.first=0.994721965713*1.00216749115;   eff.second=0.00190624905296; }
              	else if(pt<35)       { eff.first=0.993391348187*1.00229620497;   eff.second=0.00137557822021; }
              	else if(pt<40)       { eff.first=0.992284827063*1.00120405921;   eff.second=0.000934700944795; }
              	else if(pt<50)       { eff.first=0.991870039182*1.00023146204;   eff.second=0.000529497437652; }
              	else if(pt<60)       { eff.first=0.99501006222*0.999899317973;   eff.second=0.00130197167428; }
              	else if(pt<90)       { eff.first=0.990406054517*1.00095883519;   eff.second=0.00204104371381; }
              	else if(pt<140)      { eff.first=1.0090275981*1.00119553919;   eff.second=0.00635242794538; }
              	else if(pt<300)      { eff.first=1.01095605627*1.00403666587;   eff.second=0.0336631303287; }
          }
          else if(fabs(eta)<2.1)
          {
              	if(pt<20)	     { eff.first=1.01801842846*0.982435458763;  eff.second=0.00406425740549; }
              	else if(pt<25)       { eff.first=1.00035139943*0.994891457975;   eff.second=0.00172470517427; }
              	else if(pt<30)       { eff.first=0.998486035922*1.00249350572;   eff.second=0.00107768490816; }
              	else if(pt<35)       { eff.first=0.996558450006*1.00290702298;   eff.second=0.00086011870865; }
              	else if(pt<40)       { eff.first=0.996026448929*1.00192117285;   eff.second=0.000669738279806; }
              	else if(pt<50)       { eff.first=0.996061812631*1.00005805573;   eff.second=0.000274026686538; }
              	else if(pt<60)       { eff.first=0.995182727487*1.00022513914;   eff.second=0.000959995162155; }
              	else if(pt<90)       { eff.first=0.992486181067*1.00024336284;   eff.second=0.00159233268263; }
              	else if(pt<140)      { eff.first=1.02312938008*0.999888742654;   eff.second=0.00544318121633; }
              	else if(pt<300)      { eff.first=0.974754171944*0.997353820859;   eff.second=0.0296140483912; }
          }
          else if(fabs(eta)<2.4)
          {
              	if(pt<20)            { eff.first=1.00504433325*1.07619697103;  eff.second=0.008455767881; }
              	else if(pt<25)       { eff.first=0.998089082654*1.06032572307;   eff.second=0.00403060306004; }
              	else if(pt<30)       { eff.first=0.996182890451*1.04687216115;   eff.second=0.00256516752568; }
              	else if(pt<35)       { eff.first=1.00055105199*1.03180545859;   eff.second=0.0020235231465; }
              	else if(pt<40)       { eff.first=0.992563418852*1.02254502582;   eff.second=0.00166072070227; }
              	else if(pt<50)       { eff.first=0.995144128208*1.01054774978;   eff.second=0.00110937715643; }
              	else if(pt<60)       { eff.first=0.993590319967*1.00687446126;   eff.second=0.00281344868031; }
              	else if(pt<90)       { eff.first=0.9894841861*1.00450315202;   eff.second=0.00515045546188; }
              	else if(pt<140)      { eff.first=1.06017334329*0.999247787261;   eff.second=0.013319888839; }
              	else if(pt<300)      { eff.first=0.890546814737*1.00520693878;   eff.second=0.143219475222; }
          }
      
      }
	}
	break;
      }
      
      return eff;
    }

 private:

};


#endif
