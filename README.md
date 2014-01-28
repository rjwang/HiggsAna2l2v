HiggsAna2l2v
============

git clone git@github.com:rjwang/HiggsAna2l2v.git 



------------------------
scramv1 project CMSSW CMSSW_5_3_14
cd CMSSW_5_3_14/src
cmsenv

git cms-addpkg DataFormats/ParticleFlowCandidate
git cms-addpkg DataFormats/TrackReco
git cms-addpkg DataFormats/VertexReco

set FROM = "YOUR OLD RELEASE/src/"
mkdir -p CMGTools

cp -r $FROM/CMGTools/External CMGTools/
cp -r $FROM/CommonTools ./
cp -r $FROM/EGamma ./
cp -r $FROM/HiggsAnalysis ./
cp -r $FROM/pharris ./
cp -r $FROM/PhysicsTools ./
cp -r $FROM/RecoBTag ./
cp -r $FROM/RecoMET ./
cp -r $FROM/RecoParticleFlow ./
cp -r $FROM/RecoVertex ./

cd CMGTools
git clone git@github.com:rjwang/HiggsAna2l2v.git
cd -

scramv1 b -j 8 
------------------------
