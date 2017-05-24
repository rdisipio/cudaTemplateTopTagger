// Dear emacs, this is -*- c++ -*-

#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"

#include <cuda.h>
#include <math_constants.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>

using namespace std;

#define TOP_MASS 172.5
#define W_MASS    80.4
#define B_MASS     5.5
#define THREADS_PER_BLOCK 1

///////////////////////////////////////////

__global__ 
void setup_rng( curandState * state )
{
  int id = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
  curand_init( 1234, id, 0, &state[id] );
}

///////////////////////////////////////////

__host__ __device__
float Pseudorapidity( float4 p )
{
  //  if( (p.x == 0) && (p.y == 0) && (p.z == 0.) ) return 0.;
  const float P = sqrtf( p.x*p.x + p.y*p.y + p.z*p.z );
  return atanhf( p.z / P );
}


__host__ __device__
float Phi( float px, float py )
{
   return ( px == 0. && py == 0. ) ? 0. : atan2f( py, px );
}

__host__ __device__
float Phi_mpi_pi( float x )
{
   while( x >= CUDART_PI_F ) x -= 2.*CUDART_PI_F;
   while( x < -CUDART_PI_F ) x += 2.*CUDART_PI_F;

   return x;
}

__host__ __device__
float safe_pseudorapidity( float x, float abs_eta_max = 2.0 )
{
  if( x > abs_eta_max ) return abs_eta_max;
  if( x < -abs_eta_max ) return -abs_eta_max;

  return x;
}

__host__ __device__
float DeltaR_PtEtaPhiE( float4 p1, float4 p2 )
{
  const float dEta = p1.y - p2.y;
  const float dPhi = Phi_mpi_pi( p1.z - p2.z );
  float dR   = dEta*dEta + dPhi*dPhi;
  dR = ( dR > 0. ) ? sqrtf(dR) : -sqrtf(-dR);
  //  printf( "dEta=(%f-%f)=%f, dPhi=(%f-%f)=%f, dR=%f\n", p1.y, p2.y, dEta, p1.z, p2.z, dPhi, dR );
  return dR;
}

///////////////////////////////////////////


__host__ __device__
float4 PtEtaPhiE_to_PxPyPzE( float4 p )
{
  float4 v;
  v.x = p.x * cosf( p.z );
  v.y = p.x * sinf( p.z );
  v.z = p.x * sinhf( p.y );
  v.w = p.w;

  return v;
}

__host__ __device__
float4 PxPyPzE_to_PtEtaPhiE( float4 p )
{
  const float P = sqrt( p.x*p.x + p.y*p.y + p.z*p.z );

  float4 v;
  v.x = sqrt( p.x*p.x + p.y*p.y );
  v.y = 0.5 * logf( (P+p.z) / (P-p.z) );
  v.z = Phi_mpi_pi( Phi( p.x, p.y ) );
  v.w = p.w;

  return v;
}

///////////////////////////////////////////

__host__ __device__
float4 boost_PxPyPzE( const float4 p, const float4 parent )
{
  // see: https://root.cern.ch/doc/master/TLorentzVector_8cxx_source.html#l00302

  // printf( "Boost:   px=%f py=%f pz=%f E=%f\n", p.x, p.y, p.z, p.w );

  float betax = parent.x / parent.w;
  float betay = parent.y / parent.w;
  float betaz = parent.z / parent.w;
  float beta2 = betax*betax + betay*betay + betaz*betaz;
  float gamma = 1.0 / sqrtf( 1. - beta2 );
  float bp    = betax*p.x + betay*p.y + betaz*p.z;
  float gamma2 = ( beta2 > 0. ) ? ( gamma - 1. ) / beta2 : 0.; 
  //  printf( "Boost: betax=%f betay=%f betaz=%f beta2=%f gamma=%f, gamma2=%f\n", betax, betay, betaz, beta2, gamma, gamma2 );

  float4 p_prime;
  p_prime.x = p.x + gamma2*bp*betax + gamma*betax*p.w;
  p_prime.y = p.y + gamma2*bp*betay + gamma*betay*p.w;
  p_prime.z = p.z + gamma2*bp*betaz + gamma*betaz*p.w;
  p_prime.w = gamma * ( p.w + bp );

  //  printf( "Boosted: px=%f py=%f pz=%f E=%f\n", p_prime.x, p_prime.y, p_prime.z, p_prime.w );

  return p_prime;
  //  *v = host;
}

///////////////////////////////////////////

__host__ __device__
float partonOverlap_PtEtaPhiE( float4 * clusters, int n_clusters, float4 p_PtEtaPhiE, float dR_max )
{
  float ov = 0.;
  const float sigma = p_PtEtaPhiE.w / 3.;

  float E_calo = 0.;
  for( int i = 0 ; i < n_clusters ; ++i ) {
    float4 cl_PtEtaPhiE = clusters[i];

    float dR = DeltaR_PtEtaPhiE( cl_PtEtaPhiE, p_PtEtaPhiE );
    //printf("dR=%f\n", dR );
    if( dR > dR_max ) continue;

    E_calo += cl_PtEtaPhiE.w;
  }
  if( E_calo == 0. ) return 0.;

  const float dE = p_PtEtaPhiE.w - E_calo;
  ov = dE / sigma;

  return ov*ov;
}


__host__ __device__
float logOverlap_PtEtaPhiE( float4 * clusters, int n_clusters, float4 _b_PtEtaPhiE, float4 _q1_PtEtaPhiE, float4 _q2_PtEtaPhiE )
{
  float ov3 = 0.;
  float dR_max = 0.2;

  ov3 -= partonOverlap_PtEtaPhiE( clusters, n_clusters,  _b_PtEtaPhiE, dR_max );
  ov3 -= partonOverlap_PtEtaPhiE( clusters, n_clusters, _q1_PtEtaPhiE, dR_max );
  ov3 -= partonOverlap_PtEtaPhiE( clusters, n_clusters, _q2_PtEtaPhiE, dR_max );

  return ov3;
}


__device__
float4 smear_parton_PtEtaPhiE( float4 jet,  curandState * localState ) 
{
  const float j_pt  = jet.x;
  const float j_eta = jet.y;
  const float j_phi = jet.z;

  // top (eta,phi) must match geometrically the jet
  const float t_eta = safe_pseudorapidity( 0.2*curand_normal( localState ) + j_eta, 2.5 );
  const float t_phi = Phi_mpi_pi( 0.2*curand_normal( localState ) + j_phi );
  const float t_pt  = fmaxf( 0., 50.*curand_normal( localState ) + j_pt ); // trust pT more than E measurement
 
  // adjust Energy accordingly, with the constraint of TOP_MASS
  const float t_px = t_pt*cosf(t_phi);
  const float t_py = t_pt*sinf(t_phi);
  const float t_pz = t_pt*sinhf(t_eta);
  const float t_E2 = t_px*t_px + t_py*t_py + t_pz*t_pz + TOP_MASS*TOP_MASS;
  if( t_E2 < 0. ) {
    // trial resulted in an unphysical situation
    return make_float4( 0., 0., 0., 0. );
  }
  const float t_E  = sqrtf( t_E2 ); 

  //  float beta2 = t_px*t_px/t_E2 + t_py*t_py/t_E2 + t_pz*t_pz/t_E2;
  // printf("smear: top beta2 = %f\n", beta2 );

  float4 _t_PtEtaPhiE = make_float4( t_pt, t_eta, t_phi, t_E );
  
  return _t_PtEtaPhiE;
}


//////////////////////////////////////////


__global__
void tagger( float4 * jet, float4 * clusters, int * n_clusters, float * Ov3, float4 * fitted_top_PtEtaPhiE, curandState * globalState )  
{
   const int tid = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
   curandState localState = globalState[tid];

   float4 _t_PtEtaPhiE = smear_parton_PtEtaPhiE( *jet, &localState );
   if( _t_PtEtaPhiE.x == 0. ) {
     Ov3[tid] = 0.;
     fitted_top_PtEtaPhiE[tid] = _t_PtEtaPhiE;
     return;
   }

   float4 _t = PtEtaPhiE_to_PxPyPzE( _t_PtEtaPhiE );

   const float tWb_theta = acosf( 2. * curand_uniform( &localState ) - 1. );
   const float tWb_phi   = Phi_mpi_pi( 2. * CUDART_PI_F * curand_uniform( &localState ) );

   const float Wqq_theta = acosf( 2. * curand_uniform( &localState ) - 1. );
   const float Wqq_phi   = Phi_mpi_pi( 2. * CUDART_PI_F * curand_uniform( &localState ) );

   globalState[tid] = localState;

   
   // decay t->Wb
   float W_E = ( TOP_MASS*TOP_MASS + W_MASS*W_MASS - B_MASS*B_MASS ) / ( 2. * TOP_MASS );
   float B_E = ( TOP_MASS*TOP_MASS - W_MASS*W_MASS + B_MASS*B_MASS ) / ( 2. * TOP_MASS );
   float _P  = sqrtf( W_E*W_E - W_MASS*W_MASS ); 

   float _px = _P * sinf(tWb_theta) * cosf(tWb_phi);
   float _py = _P * sinf(tWb_theta) * sinf(tWb_phi);
   float _pz = _P * cosf(tWb_theta);

   float4 _W = make_float4(  _px,  _py,  _pz, W_E );
   float4 _b = make_float4( -_px, -_py, -_pz, B_E );

   _W = boost_PxPyPzE( _W, _t );
   _b = boost_PxPyPzE( _b, _t );

   // decay W->qq
   // assume both q massless?
   float q1_E  = 0.5 * W_MASS;
   float q2_E  = 0.5 * W_MASS;
   float Wqq_P = sqrtf( q1_E*q1_E ); // -LQ_MASS*LQ_MASS
   _px = Wqq_P * sinf(Wqq_theta) * cosf(Wqq_phi);
   _py = Wqq_P * sinf(Wqq_theta) * sinf(Wqq_phi);
   _pz = Wqq_P * cosf(Wqq_theta);
   float4 _q1 = make_float4(  _px,  _py,  _pz, q1_E );
   float4 _q2 = make_float4( -_px, -_py, -_pz, q2_E );

   //printf("_q1 pre-boost: px=%f py=%f pz=%f\n", _q1.x, _q1.y, _q1.z );
   _q1 = boost_PxPyPzE( _q1, _W );
   //printf("_q1 pst-boost: px=%f py=%f pz=%f\n", _q1.x, _q1.y, _q1.z );
   _q2 = boost_PxPyPzE( _q2, _W );


   float4 _b_PtEtaPhiE  = PxPyPzE_to_PtEtaPhiE( _b );
   float4 _q1_PtEtaPhiE = PxPyPzE_to_PtEtaPhiE( _q1 );
   float4 _q2_PtEtaPhiE = PxPyPzE_to_PtEtaPhiE( _q2 );

   float logOv3 = logOverlap_PtEtaPhiE( clusters, *n_clusters, _b_PtEtaPhiE, _q1_PtEtaPhiE, _q2_PtEtaPhiE );

   Ov3[tid] = logOv3;
   fitted_top_PtEtaPhiE[tid] = _t_PtEtaPhiE;   
}

///////////////////////////////////////////


void run_tagger( int n_points, float4 * jet, float4 * clusters, int * n_clusters, float * Ov3, float4 * fitted_top_PtEtaPhiE, curandState * globalState )
{
  const int n_blocks = n_points / THREADS_PER_BLOCK;

  tagger<<< n_blocks, THREADS_PER_BLOCK >>>( jet, clusters, n_clusters, Ov3, fitted_top_PtEtaPhiE, globalState );

  cudaDeviceSynchronize();
}


///////////////////////////////////////////


int main( int argc, char ** argv )
{ 
   int success = 0;

   gROOT->SetBatch();

   const char * infilename = argv[1];

//   Int_t             eventNumber;
   vector< float > * jet_pt = NULL;
   vector< float > * jet_eta = NULL;
   vector< float > * jet_phi = NULL;
   vector< float > * jet_E = NULL;
   vector< float > * jet_m = NULL;
   vector< int >   * cl_matched_jet = NULL;
   vector< float > * cl_pt = NULL;
   vector< float > * cl_eta = NULL;
   vector< float > * cl_phi = NULL;
   vector< float > * cl_E = NULL;

   TBranch * b_jet_pt = NULL;
   TBranch * b_jet_eta = NULL;
   TBranch * b_jet_phi = NULL;
   TBranch * b_jet_E = NULL;
   TBranch * b_jet_m = NULL;
   TBranch * b_cl_matched_jet = NULL;
   TBranch * b_cl_pt  = NULL;
   TBranch * b_cl_eta = NULL;
   TBranch * b_cl_phi = NULL;
   TBranch * b_cl_E   = NULL;


   TChain tree( "physics" );
   tree.Add( infilename );

 //  tree.SetBranchAddress( "eventNumber", &eventNumber );
   tree.SetBranchAddress( "jet_pt",      &jet_pt, &b_jet_pt );
   tree.SetBranchAddress( "jet_eta",     &jet_eta, &b_jet_eta );
   tree.SetBranchAddress( "jet_phi",     &jet_phi, &b_jet_phi );
   tree.SetBranchAddress( "jet_E",       &jet_E, &b_jet_E );
   tree.SetBranchAddress( "jet_m",       &jet_m, &b_jet_m );
   tree.SetBranchAddress( "cl_matched_jet", &cl_matched_jet, &b_cl_matched_jet );
   tree.SetBranchAddress( "cl_pt",   &cl_pt, &b_cl_pt );
   tree.SetBranchAddress( "cl_eta",  &cl_eta, &b_cl_eta );
   tree.SetBranchAddress( "cl_phi",  &cl_phi, &b_cl_phi );
   tree.SetBranchAddress( "cl_E",    &cl_E, &b_cl_E );

   TFile * f_outfile = TFile::Open( "Ov3.root", "RECREATE" );
   TH1F * h_ov3 = new TH1F( "Ov3", "Overlap", 20, 0., 1. );

   // GPU Memory allocation
   float4 * gpu_jet_fourmom = NULL;
   float4 * gpu_fit_top     = NULL;
   float4 * gpu_cl_fourmom  = NULL;
   int    * gpu_n_clusters  = NULL;
   float  * gpu_Ov3         = NULL;

   // GPU random numbers generator
   const int n_points = 1024;
   const int n_blocks = n_points / THREADS_PER_BLOCK;
   curandState  * gpu_rng_states = new curandState[n_points];
   cudaMallocManaged( &gpu_rng_states, n_points*sizeof(curandState) );

   cudaMallocManaged( &gpu_jet_fourmom, sizeof(float4) );
   cudaMallocManaged( &gpu_Ov3, n_points * sizeof(float) );
   cudaMallocManaged( &gpu_fit_top, n_points*sizeof(float4) );
   cudaMallocManaged( &gpu_n_clusters, sizeof(int) );
 
   setup_rng<<< n_blocks, THREADS_PER_BLOCK >>>( gpu_rng_states );

   cudaDeviceSynchronize();
   cout << "INFO: RNG initialization done" << endl;


   // Start reading event tree
   const Int_t nEntries = tree.GetEntries();
   cout << "INFO: number of entries: " << nEntries << endl;

   for( Int_t ientry = 0 ; ientry < nEntries ; ++ientry ) {
     Long64_t centry = tree.GetEntry( ientry );

     const Int_t n_jets     = jet_pt->size();
     const Int_t n_clusters = cl_pt->size();

     cout << "INFO: event: " << ientry << " N jets = " << n_jets << " N clusters = " << n_clusters << endl;

     for( Int_t ijet = 0 ; ijet < n_jets ; ++ijet ) {
       if( jet_pt->at(ijet) == 0 ) continue;

       gpu_jet_fourmom[0].x = jet_pt->at(ijet);
       gpu_jet_fourmom[0].y = jet_eta->at(ijet);
       gpu_jet_fourmom[0].z = jet_phi->at(ijet);
       gpu_jet_fourmom[0].w = jet_E->at(ijet);

       cout << "Jet pt=" << gpu_jet_fourmom[0].x << " eta=" << gpu_jet_fourmom[0].y << " phi=" 
	    << gpu_jet_fourmom[0].z << " E=" << gpu_jet_fourmom[0].w << endl;

       // count matched clusters
       int n_cl_jet = 0;
       for( int k = 0 ; k < n_clusters ; ++k ) {
         if( cl_matched_jet->at(k) != ijet ) continue;
         n_cl_jet++;
       }

       // reserve memory on GPU
       cudaMallocManaged((void**)&gpu_cl_fourmom,  n_cl_jet*sizeof(float4) );
       gpu_n_clusters[0] = n_cl_jet;

       int icluster = 0;
       for( Int_t k = 0 ; k < n_clusters ; ++k ) {
          if( cl_matched_jet->at(k) != ijet ) continue;

          // fill clusters vector
          float4 v;
          v.x = cl_pt->at(k);
          v.y = cl_eta->at(k);
          v.z = cl_phi->at(k);
          v.w = cl_E->at(k);       
          gpu_cl_fourmom[icluster] = v;
	  icluster++;
       } // clusters loop

       run_tagger( n_points, gpu_jet_fourmom, gpu_cl_fourmom, gpu_n_clusters, gpu_Ov3, gpu_fit_top, gpu_rng_states );

       // retrieve result
       float  bestfit_ov3 = -1e11;
       float4 bestfit_top;
       for( int ipoint = 0 ; ipoint < n_points ; ++ipoint ) {
	 float ov3 = gpu_Ov3[ipoint];
	 if( ov3 == 0. ) continue;

	 if( ov3 < bestfit_ov3 ) continue;

	 bestfit_ov3 = ov3;
	 bestfit_top = gpu_fit_top[ipoint];
	 
       }
       cout << "Best-fit Ov3=" << bestfit_ov3 << " : t_pt=" << bestfit_top.x << " eta=" << bestfit_top.y 
	    << " phi=" << bestfit_top.z << " E=" << bestfit_top.w << endl;

       if( bestfit_ov3 != -1e11 ) {
	 h_ov3->Fill( exp(bestfit_ov3) );
       }
       else {
	 h_ov3->Fill( -1 );
       }

       cudaFree( gpu_cl_fourmom );

     } // jets loop


   } // event loop

   cudaFree( gpu_Ov3 );
   cudaFree( gpu_rng_states );

   f_outfile->cd();
   h_ov3->Write();
   f_outfile->Close();

   return success;
}
