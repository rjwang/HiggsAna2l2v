#ifndef MuonTriggerEfficiencySF_h
#define MuonTriggerEfficiencySF_h

/** \class MuonTriggerEfficiencySF_h
 *  No description available.
 *
 *  $Date: 2014/03/14 14:11:23 $
 *  $Revision: 1.10 $
 *  \author Renjie Wang
 */

class MuonTriggerEfficiencySF {
public:

    //
    MuonTriggerEfficiencySF() { }

    //
    ~MuonTriggerEfficiencySF() {}

    //
    std::pair<float,float> getMuonTriggerEfficiencySF(float eta1, float eta2) {

        std::pair<float,float> eff(1.0,0.0);
        if(eta1 > -2.4 && eta1 < -1.2) {

            if     (eta2 >  -2.4 && eta2 < -1.2)     { eff.first =1.0067 ; eff.second =0.0104 ; } 
	    else if(eta2 >= -1.2 && eta2 < 0   )     { eff.first =0.9983 ; eff.second =0.0108 ; } 
	    else if(eta2 >= 0    && eta2 < 1.2 )     { eff.first =1.0164 ; eff.second =0.0110 ; } 
	    else if(eta2 >= 1.2  && eta2 < 2.4 )     { eff.first =1.0868 ; eff.second =0.0112 ; }

        } else if(eta1 >= -1.2 && eta1 < 0) {

            if     (eta2 >  -2.4 && eta2 < -1.2)     { eff.first =0.9833 ; eff.second =0.0106 ; } 
	    else if(eta2 >= -1.2 && eta2 < 0   )     { eff.first =0.9885 ; eff.second =0.0112 ; } 
	    else if(eta2 >= 0    && eta2 < 1.2 )     { eff.first =0.9870 ; eff.second =0.0112 ; } 
	    else if(eta2 >= 1.2  && eta2 < 2.4 )     { eff.first =0.9908 ; eff.second =0.0108 ; }

        } else if(eta1 >= 0 && eta1 < 1.2) {

            if     (eta2 >  -2.4 && eta2 < -1.2)     { eff.first =1.0015 ; eff.second =0.0108 ; } 
	    else if(eta2 >= -1.2 && eta2 < 0   )     { eff.first =0.9874 ; eff.second =0.0112 ; } 
	    else if(eta2 >= 0    && eta2 < 1.2 )     { eff.first =0.9915 ; eff.second =0.0113 ; } 
	    else if(eta2 >= 1.2  && eta2 < 2.4 )     { eff.first =0.9964 ; eff.second =0.0108 ; }

        } else if(eta1 >= 1.2 && eta1 < 2.4) {

            if     (eta2 >  -2.4 && eta2 < -1.2)     { eff.first =1.0869 ; eff.second =0.0112 ; } 
	    else if(eta2 >= -1.2 && eta2 < 0   )     { eff.first =1.0007 ; eff.second =0.0109 ; } 
	    else if(eta2 >= 0    && eta2 < 1.2 )     { eff.first =0.9998 ; eff.second =0.0109 ; } 
	    else if(eta2 >= 1.2  && eta2 < 2.4 )     { eff.first =1.0075 ; eff.second =0.0105 ; }

        }

        return eff;

    } //getMuonTriggerEfficiencySF

private:

};


#endif
