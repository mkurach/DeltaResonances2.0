/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
 *                                                                              *
 * Version:                                                                     *
 *      Release, 2.0.3, 1 February 2011                                         *
 *                                                                              *
 * Authors:                                                                     *
 *      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
 *      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
 *      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
 *      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
 *                                                                              *
 * Project homepage:                                                            *
 *      http://therminator2.ifj.edu.pl/                                         *
 *                                                                              *
 * For the detailed description of the program and further references           *
 * to the description of the model please refer to                              *
 * http://arxiv.org/abs/1102.0273                                               *
 *                                                                              *
 * This code can be freely used and redistributed. However if you decide to     *
 * make modifications to the code, please, inform the authors.                  *
 * Any publication of results obtained using this code must include the         *
 * reference to arXiv:1102.0273 and the published version of it, when           *
 * available.                                                                   *
 *                                                                              *
 ********************************************************************************/

#include <sys/stat.h>
#include "Configurator.h"
#include "Model.h"
#include "Crc32.h"
#include "THGlobal.h"
#include <map>
#include <vector>

using namespace std;
using namespace TMath;

extern Configurator *sMainConfig;
extern TString	sEventDIR;
extern int	sModel;
extern int	sIntegrateSample;
extern int	sRandomize;

Model::Model()
    : Xt(0.0), Xx(0.0), Xy(0.0), Xz(0.0),
    Pe(0.0), Px(0.0), Py(0.0), Pz(0.0),
    mHyperCube(0.0),
    mRandom(0),
    mSpectralFncs(0), mSpectralFncIntegrals(0),
    mDB(0)
{
  mName="";
  mHash="";
  mDescription="";
}

Model::Model(TRandom2* aRandom)
: Xt(0.0), Xx(0.0), Xy(0.0), Xz(0.0),
Pe(0.0), Px(0.0), Py(0.0), Pz(0.0),
mHyperCube(0.0),
mRandom(aRandom)
{
  mName="";
  mHash="";
  mDescription="";
  mSpectralFncs = new std::map<int, TF1*>;
  mSpectralFncIntegrals = new std::map<int, double>;

}
Model::~Model()
{
  delete mSpectralFncs;
  delete mSpectralFncIntegrals;
}

void Model::AddParameterBranch(TTree* aTree)
{
  Model_t tPar;
  
  tPar.dummy = -1.0;
  aTree->Branch(_MODEL_T_BRANCH_, &tPar, _MODEL_T_FORMAT_);
}

void Model::SetParticlePX(Particle* aParticle)
{
  aParticle->SetParticlePX(Pe,Px,Py,Pz,Xt,Xx,Xy,Xz);
}  

double Model::GetHyperCubeVolume()
{
  return mHyperCube;
}

const char* Model::GetHash()
{
  return mHash.Data();
}

const char* Model::GetName()
{
  return mName.Data();
}

const char* Model::GetDescription()
{
  return mDescription.Data();
}

void Model::CreateEventSubDir()
{
  struct stat tStatus;
  TString tEventDir = sMainConfig->GetParameter("EventDir"); tEventDir.Prepend("./");
  TString tSubDirs  = sEventDIR;
  int     tPos      = 0;

  tSubDirs.ReplaceAll(tEventDir,"");
  while(tPos <= tSubDirs.Length()) {   
    if(tEventDir.EndsWith("/")) {
      if (stat(tEventDir.Data(), &tStatus) == -1) {
	PRINT_DEBUG_3("<Model::CreateEventSubDir>\tDirectory " <<tEventDir<<" does not exist. Tying to create.");
        if (mkdir(tEventDir.Data(), S_IRWXU | S_IXUSR | S_IRWXG | S_IXGRP | S_IROTH | S_IXOTH) == -1) {
          PRINT_MESSAGE("<Model::CreateEventSubDir>\tFailed to create directory " << tEventDir);
          exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
        } else {
          PRINT_DEBUG_2("<Model::CreateEventSubDir>\tDirectory "<<tEventDir<<" created.");
        }
      } else {
	PRINT_DEBUG_3("<Model::CreateEventSubDir>\tDirectory " <<tEventDir<<" does exist.");
      }
    }
    tEventDir += tSubDirs[tPos];
    tPos++;
  }  
}

void Model::CalculateHash(TString aPreHash) {
  Crc32 tHash;
  tHash.Update(aPreHash.Data(), aPreHash.Length());
  tHash.Finish();
  mHash = tHash.GetValueHex();
  PRINT_DEBUG_1("<Model::Hash>\t"<<aPreHash.Data()<<" -> 0x"<<mHash);
}


TH1F * p33() {

    static bool create = true;

    // this is from PML eBW (from his data33 column marked as A)
    static Double_t y4[] = {

	0.000335537, 0.00200371, 0.00464069, 0.00821601,
	0.0127918, 0.0184775, 0.0254215, 0.0338109, 0.0438759,
	0.0558947, 0.0702011, 0.0871924, 0.107337, 0.13118,
	0.159347, 0.192535, 0.231494, 0.276977, 0.329664,
	0.390015, 0.458079, 0.533225, 0.613851, 0.697139,
	0.778991, 0.854333, 0.917805, 0.964779, 0.992344,
	0.999881, 0.988991, 0.962874, 0.925482, 0.880741,
	0.832059, 0.782097, 0.732765, 0.685319, 0.640506,
	0.598705, 0.56004, 0.524474, 0.49187, 0.462039,
	0.434765, 0.40983, 0.387017, 0.366124, 0.346964,
	0.329365, 0.313172, 0.298246, 0.284462, 0.271709,
	0.259888, 0.24891, 0.238697, 0.229178, 0.220292,
	0.21198, 0.204195, 0.19689, 0.190026, 0.183567,
	0.177479, 0.171733, 0.166303, 0.161166, 0.156298,
	0.151681, 0.147296, 0.143128, 0.13916, 0.13538,
	0.131776, 0.128334, 0.125046, 0.121902, 0.118892,
	0.116009, 0.113244, 0.110591, 0.108044, 0.105596,
	0.103243, 0.100977, 0.0987961, 0.0966944, 0.094668,
	0.092713, 0.0908259, 0.0890032, 0.0872417, 0.0855385,
	0.0838908, 0.0822961, 0.0807517, 0.0792555, 0.0778053,
	0.076399, 0.0750347, 0.0737106, 0.0724249, 0.0711761,
	0.0699626, 0.0687831, 0.067636, 0.0665201, 0.0654342,
	0.0643771, 0.0633477, 0.062345, 0.0613679, 0.0604156,
	0.059487, 0.0585814, 0.0576979, 0.0568357, 0.0559942,
	0.0551725, 0.05437, 0.0535861, 0.0528202, 0.0520716,
	0.0513397, 0.0506241, 0.0499242, 0.0492396, 0.0485696,
	0.047914, 0.0472722, 0.0466439, 0.0460285, 0.0454259
    } ;
    // original from PML
    static Double_t y5[] = {
        0.487868830660434  ,  0.974036309173992  ,   1.48796409314226  ,   1.88297294151952  ,   2.28387966214487  ,   2.69710433117348  ,   3.06126905224752  ,   3.48369820687455  ,   3.94681509541577  ,   4.40057706566204  ,
        4.93936605645115  ,   5.48283889953534  ,   6.11314163074160  ,   6.81950751171442  ,   7.55858848427632  ,   8.39291146773919  ,   9.34237968621590  ,   10.3691236214243  ,   11.4890241422109  ,   12.6864224931698  ,
        13.9961542566581  ,   15.3159690791888  ,   16.5829494782608  ,   17.8292131736408  ,   18.9152955911726  ,   19.7159202944624  ,   20.2343075893699  ,   20.3069956416159  ,   20.0125696007672  ,   19.4358933652416  ,
        18.4928474508350  ,   17.3823315861657  ,   16.1716157278873  ,   14.9038674730610  ,   13.6447991751455  ,   12.5056812406078  ,   11.4205163366001  ,   10.4059676207370  ,   9.56824474209078  ,   8.71912443140713  ,
        8.01315763782008  ,   7.39793576023726  ,   6.80293967011530  ,   6.32923543881047  ,   5.89510052808095  ,   5.48821214622965  ,   5.13744465308651  ,   4.85174031833841  ,   4.57035038945562  ,   4.30429423098407  ,
        4.07695313297681  ,   3.89130969244800  ,   3.71215490911060  ,   3.56308036882298  ,   3.40112470196175  ,   3.25304206818625  ,   3.16482535139844  ,   3.04715152970495  ,   2.96881716587616  ,   2.88140869791066  ,
        2.78585022621101  ,   2.72700023847186  ,   2.64608188068563  ,   2.60279320128658  ,   2.54784690812262  ,   2.49439864564337  ,   2.44898753815908  ,   2.37474004226641  ,   2.38889146786364  ,   2.32064636249309  ,
        2.25598151655673  ,   2.26163499621386  ,   2.21317590725296  ,   2.18258711838417  ,   2.14657132226298  ,   2.11184984790125  ,   2.07710318841325  ,   2.04127250297028  ,   2.00975961369809  ,   1.96467048836692  ,
        1.96033089737522  ,   1.92523920165693  ,   1.86071345511037  ,   1.83261799374429  ,   1.80431223657368  ,   1.74073722854448  ,   1.70206527339233  ,   1.70885460224525  ,   1.62966568770821  ,   1.60909660760260  ,
        1.57063349526692  ,   1.53328104700055  ,   1.50929453898890  ,   1.44028961773135  ,   1.41651195253619  ,   1.37802510806219  ,   1.34406780855180  ,   1.30897465984558  ,   1.27409180694395  ,   1.23918522190402  ,
        1.20425635771394  ,   1.16951260420208  ,   1.13384329729916  ,   1.10156913915201  ,   1.05909016156322  ,   1.04733955656242  ,   1.03466485115872  ,  0.994816847474634  ,  0.955153954468695  ,  0.941716527201837  ,
        0.932576069943827  ,  0.882687531334134  ,  0.872806631363069  ,  0.862001630988857  ,  0.814744066283836  ,  0.798214874167238  ,  0.807009326536368  ,  0.792949729784183  ,  0.738143504668653  ,  0.750470752195434  ,
        0.762983110400590  ,  0.707412710433718  ,  0.695986993562356  ,  0.697185321244401  ,  0.703998672833197  ,  0.678025155222950  ,  0.665835263500309  ,  0.640255699384825  ,  0.644986822224574  ,  0.654246037039566  ,
        0.617376174539702  ,  0.639099508356292  ,  0.603315124799318  ,  0.608256543443586  , 
    };

    static TH1F *h5 = new TH1F("hBform","B(m)",134,1.078,1.743); // PML paper
    static TH1F *h4 = new TH1F("hV","V(m)",134,1.073,1.748); // vova

    if (create) {
        for (int b = 1; b <= h5->GetNbinsX(); ++b) {
            h5->SetBinContent(b,y5[b-1]);
            h4->SetBinContent(b,y4[b-1]);
        }
        create = false;
    }

    //    return h4;  // eBW
    return h5;   // PML

}


TH1F *getBreitWigner(int pdg) {
	static std::map<int, std::vector<Double_t>> histograms;
	static std::map<int, Double_t> mMin;
	static std::map<int, Double_t> mMax;
	histograms[12112] = {0.326505,0.3539,0.384105,0.417456,0.454328,0.495138,0.540348,0.590459,0.646012,0.707571,0.775712,0.850986,0.933876,1.02473,1.12367,1.23045,1.34435,1.46395,1.58699,1.71026,1.82954,1.93977,2.03543,2.11109,2.16218,2.18565,2.18054,2.14805,2.09129,2.01471,1.9234,1.82242,1.71634,1.60895,1.50317,1.40109,1.3041,1.213,1.12817,1.04967,0.977346,0.910923,0.850037,0.794286,0.74326,0.696554,0.653783,0.614583,0.578622,0.545592,0.515215,0.487239,0.461437,0.437605,0.415558,0.395133,0.376181,0.35857,0.342182,0.326908,0.312653,0.299332,0.286864,0.275181,0.26422,0.253921,0.244235,0.235113,0.226512,0.218395,0.210725,0.20347,0.196601,0.190091,0.183915,0.17805,0.172476,0.167174,0.162125,0.157315,0.152728,0.14835,0.144169,0.140173,0.136351,0.132692,0.129188,0.125829,0.122608,0.119517,0.11655,0.113698,0.110957,0.10832,0.105782,0.103339,0.100985,0.0987155,0.0965272,0.0944159};
	mMin[12112] = 1.07454;
	mMax[12112] = 2.49;
	histograms[1214] = {0.0740463,0.0774495,0.0810549,0.0848782,0.0889366,0.0932492,0.0978367,0.102722,0.107931,0.113492,0.119435,0.125796,0.132612,0.139928,0.147791,0.156256,0.165382,0.175239,0.185904,0.197463,0.210017,0.223678,0.238576,0.254859,0.272697,0.292286,0.313855,0.337666,0.364028,0.393302,0.425913,0.462359,0.503236,0.549247,0.601238,0.660218,0.727406,0.804273,0.892601,0.994551,1.11275,1.25036,1.41121,1.59985,1.82155,2.08228,2.38825,2.74517,3.15658,3.62099,4.12769,4.65154,5.14911,5.56062,5.82197,5.88621,5.74336,5.42462,4.98786,4.49539,3.99766,3.52747,3.10178,2.72641,2.40055,2.12007,1.87954,1.67342,1.49652,1.34429,1.21283,1.09884,0.999591,0.912794,0.83656,0.769319,0.709766,0.656811,0.609544,0.567198,0.52913,0.494793,0.463723,0.435525,0.40986,0.386437,0.365004,0.345343,0.327266,0.310608,0.295224,0.280988,0.267789,0.255528,0.244119,0.233484,0.223554,0.214268,0.205572,0.197416};
	mMin[1214] = 1.07454;
	mMax[1214] = 1.88;
	histograms[22112] = {0.0860304,0.0903747,0.0950003,0.0999312,0.105194,0.110817,0.116832,0.123277,0.130191,0.137618,0.145608,0.154218,0.163511,0.173558,0.184438,0.196242,0.209075,0.223053,0.23831,0.255001,0.273301,0.293415,0.315578,0.340063,0.367188,0.397325,0.430908,0.46845,0.510554,0.557935,0.611441,0.672081,0.741055,0.819796,0.910011,1.01373,1.13333,1.27163,1.43182,1.61747,1.83236,2.08011,2.36349,2.68329,3.0365,3.41394,3.79769,4.15969,4.46332,4.66993,4.74949,4.69096,4.50639,4.22588,3.88708,3.52505,3.16651,2.82861,2.52036,2.24502,2.00228,1.78992,1.60486,1.44378,1.30351,1.18115,1.07415,0.980297,0.897701,0.824761,0.760121,0.702634,0.651331,0.605391,0.564118,0.52692,0.493293,0.462803,0.435081,0.409808,0.386707,0.36554,0.346099,0.328204,0.311695,0.296435,0.2823,0.269184,0.256991,0.245635,0.235043,0.225148,0.215888,0.207211,0.199068,0.191417,0.184218,0.177436,0.17104,0.164999};
	mMin[22112] = 1.07454;
	mMax[22112] = 1.985;
	histograms[32214] = {0.146724,0.157887,0.170082,0.183426,0.198056,0.214126,0.231813,0.251317,0.272868,0.29673,0.323202,0.352629,0.385401,0.421964,0.46282,0.508535,0.559735,0.617105,0.681377,0.753303,0.833613,0.922945,1.02173,1.13006,1.24744,1.37258,1.50306,1.63512,1.76356,1.88185,1.98268,2.05891,2.10475,2.11689,2.09516,2.0425,1.96427,1.8671,1.75779,1.64251,1.52623,1.41267,1.30432,1.20269,1.10854,1.02206,0.943125,0.871354,0.80626,0.747296,0.693905,0.645549,0.601718,0.561943,0.525798,0.4929,0.462905,0.435508,0.410436,0.38745,0.366336,0.346906,0.328991,0.312444,0.297133,0.282942,0.269765,0.257511,0.246097,0.235448,0.225499,0.21619,0.207468,0.199284,0.191595,0.184362,0.17755,0.171127,0.165063,0.159332,0.153909,0.148774,0.143905,0.139285,0.134897,0.130724,0.126754,0.122972,0.119368,0.115929,0.112647,0.10951,0.106511,0.103642,0.100894,0.0982618,0.0957378,0.0933164,0.0909919,0.088759};
	mMin[32214] = 1.07454;
	mMax[32214] = 2.65;
	histograms[2122] = {0.057436,0.0602346,0.0632029,0.0663542,0.0697033,0.0732663,0.077061,0.0811071,0.0854263,0.0900426,0.094983,0.100277,0.105958,0.112063,0.118633,0.125716,0.133363,0.141634,0.150596,0.160324,0.170905,0.182438,0.195035,0.208825,0.223959,0.240607,0.25897,0.27928,0.30181,0.326877,0.354857,0.386194,0.421415,0.461148,0.506145,0.557307,0.615722,0.682706,0.759849,0.849079,0.952729,1.07361,1.21507,1.38104,1.57601,1.80482,2.0722,2.38178,2.73433,3.12492,3.53917,3.94957,4.31456,4.58408,4.71331,4.6797,4.49275,4.18921,3.81762,3.42279,3.03733,2.68055,2.36123,2.08127,1.83872,1.62988,1.45047,1.2963,1.16355,1.04892,0.949567,0.86311,0.787561,0.721264,0.662842,0.611147,0.565224,0.52427,0.487614,0.454689,0.425015,0.398186,0.373855,0.351725,0.331542,0.313084,0.296164,0.280615,0.266293,0.253074,0.240846,0.229514,0.218992,0.209204,0.200083,0.19157,0.183611,0.17616,0.169173,0.162613};
	mMin[2122] = 1.07454;
	mMax[2122] = 2.07;
	histograms[32212] = {0.0504494,0.0528862,0.0554679,0.0582056,0.0611116,0.0641994,0.0674837,0.0709809,0.0747089,0.0786876,0.0829391,0.0874878,0.0923609,0.0975888,0.103205,0.109248,0.11576,0.122789,0.130389,0.138621,0.147554,0.157267,0.16785,0.179404,0.192049,0.20592,0.221173,0.23799,0.256584,0.277201,0.300132,0.325718,0.354364,0.386549,0.422848,0.463946,0.510671,0.564018,0.625199,0.695684,0.777263,0.872117,0.982901,1.11282,1.2657,1.44604,1.65887,1.90952,2.20282,2.54166,2.92442,3.34118,3.76914,4.1696,4.4906,4.67894,4.69954,4.55108,4.26615,3.89599,3.49158,3.09156,2.71946,2.38632,2.09491,1.84339,1.62776,1.44338,1.28567,1.15048,1.03423,0.933867,0.846853,0.771076,0.704788,0.646542,0.595142,0.549593,0.509066,0.472868,0.440418,0.411226,0.384877,0.361019,0.339352,0.319617,0.301593,0.28509,0.269941,0.256004,0.243152,0.231275,0.220278,0.210076,0.200592,0.191762,0.183526,0.175832,0.168633,0.161887};
	mMin[32212] = 1.07454;
	mMax[32212] = 2.1;
	histograms[2216] = {0.0454789,0.0476624,0.0499738,0.0524228,0.0550199,0.0577769,0.0607065,0.0638227,0.0671412,0.0706789,0.0744548,0.0784898,0.0828072,0.087433,0.0923959,0.0977281,0.103466,0.10965,0.116325,0.123543,0.131362,0.139849,0.149077,0.159134,0.170116,0.182136,0.195324,0.20983,0.225827,0.243519,0.263142,0.284974,0.309343,0.336639,0.367322,0.401945,0.441172,0.485801,0.536799,0.595344,0.662869,0.741128,0.83227,0.938919,1.06428,1.2122,1.38726,1.5947,1.8402,2.12917,2.46543,2.84856,3.27002,3.70798,4.12347,4.46205,4.66603,4.69573,4.54776,4.25625,3.87609,3.4617,3.05382,2.67661,2.34088,2.04884,1.79806,1.58405,1.40178,1.24642,1.11366,0.999795,0.901727,0.816876,0.743115,0.678694,0.622169,0.572351,0.528252,0.489056,0.454079,0.422748,0.394585,0.369181,0.346193,0.325327,0.306332,0.288992,0.273121,0.25856,0.245167,0.232821,0.221417,0.210859,0.201067,0.191968,0.183497,0.175599,0.168221,0.161319};
	mMin[2216] = 1.07454;
	mMax[2216] = 2.125;
	histograms[12216] = {0.0384757,0.040217,0.0420544,0.0439947,0.0460453,0.0482145,0.050511,0.0529447,0.0555261,0.0582671,0.0611805,0.0642804,0.0675826,0.0711042,0.0748644,0.0788844,0.0831879,0.087801,0.0927532,0.0980771,0.10381,0.109992,0.116672,0.1239,0.131738,0.140254,0.149523,0.159636,0.170694,0.182815,0.196132,0.210805,0.227015,0.244977,0.264943,0.287209,0.312127,0.340116,0.371679,0.40742,0.448069,0.494516,0.547849,0.609406,0.680837,0.764192,0.862023,0.977509,1.11461,1.27824,1.47438,1.7102,1.99382,2.33365,2.73654,3.20418,3.7267,4.27404,4.7886,5.18829,5.38953,5.34528,5.07109,4.63517,4.12248,3.60298,3.11982,2.69258,2.3256,2.01522,1.75454,1.53597,1.3524,1.19767,1.06664,0.955077,0.859564,0.777326,0.706125,0.644147,0.589918,0.542234,0.500109,0.462729,0.429419,0.39962,0.372859,0.348743,0.326937,0.307156,0.28916,0.272741,0.257719,0.243941,0.231274,0.219601,0.20882,0.198842,0.18959,0.180994};
	mMin[12216] = 1.07454;
	mMax[12216] = 2.07;
	static TH1F *hist = new TH1F("hist","hist",100,mMin[pdg],mMax[pdg]);
	for(int i = 1; i < hist->GetNbinsX(); i++)
		hist->SetBinContent(i,histograms[pdg][i-1]);

	return hist;
}


double Model::CalcMass(ParticleType *aPartType, double &statWeight) {
    double m0 = aPartType->GetMass();
    double Gamma = aPartType->GetGamma();

    if (Gamma < 1e-6) {
        statWeight = 1.0;
        return m0;
    }

    double m_min = m0 - 2*Gamma;
    double m_max = m0 + 2*Gamma;
    if (m_min < 0.) {
        m_min = 0;
        m_max = 2*m0;
    }

    double norm = 1.0;  
    TF1 *fnc;
    if (mSpectralFncs->find(aPartType->GetPDGCode()) == mSpectralFncs->end()) {
        fnc = new TF1(Form("spectralFunction_pdg%i",aPartType->GetPDGCode()),"x*[0]*[1]/((x*x-[0]*[0])*(x*x-[0]*[0])+[0]*[0]*[1]*[1])",m_min,m_max);
        fnc->SetParameters(m0,Gamma);
        norm = fnc->Integral(m_min,m_max);
        
        (*mSpectralFncs)[aPartType->GetPDGCode()] = fnc;
        (*mSpectralFncIntegrals)[aPartType->GetPDGCode()] = norm;

       // cout << "mass, Gamma, m_min, m_max " << m0 << " " << Gamma << " " << m_min << " " << m_max << endl;
       // cout << "norm " << norm << endl;
    } else {
        fnc = (*mSpectralFncs)[aPartType->GetPDGCode()];
        norm = (*mSpectralFncIntegrals)[aPartType->GetPDGCode()];
    }
    double m = fnc->GetRandom(m_min,m_max);
    statWeight = fnc->Eval(m) / fnc->Eval(m0) / (m_max-m_min);
    return m;
}


double Model::CalcMassDistr(ParticleType *aPartType, double &statWeight) {

    int pdg = aPartType->GetPDGCode();
    if (Abs(pdg) == 2224 || Abs(pdg) == 2214 || Abs(pdg) == 2114 || Abs(pdg) == 1114) {
        TH1F *distr = p33();
        double m = distr->GetRandom();
        statWeight = distr->Interpolate(m);
        return m;
    } 
    
    else if(  Abs(pdg) == 12112||Abs(pdg) == 1214||Abs(pdg) == 22112||Abs(pdg) == 32214||Abs(pdg) == 2122||Abs(pdg) == 32212||Abs(pdg) == 2216||Abs(pdg) == 12216) {
      TH1F *distr = getBreitWigner(Abs(pdg));
        double m = distr->GetRandom();
        statWeight = 1;
        return m;
    }
    else {
	      statWeight = 1;
        return aPartType->GetMass();
       // return CalcMass(aPartType, statWeight);
    }


}

void Model::GetParticleMass(ParticleType *aPartType, bool finiteWidth, double &M, double &spectralFunctionWeight ) {
    if (finiteWidth) {
	M = CalcMassDistr(aPartType, spectralFunctionWeight);
    } else {
	M = aPartType->GetMass();
	spectralFunctionWeight = 1.0;
    }
}


