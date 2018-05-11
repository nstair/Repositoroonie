/*
 * @BEGIN LICENSE
 *
 * so_mp2 by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"

#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/process.h"


// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

// INSERT_ELEMENTS (collection, first, last)
// - fill values from first to last into the collection
// - NOTE: NO half-open range
template <typename T>
inline void INSERT_ELEMENTS(T &coll, int first, int last) {
  for (int i = first; i < last; ++i) {
    coll.insert(coll.end(), i);
  }
}

// PRINT_ELEMENTS()
// - prints optional string optcstr followed by
// - all elements of the collection coll
// - separated by spaces
template <typename T>
inline void PRINT_ELEMENTS(const T &coll, const std::string &optcstr = "",
                           bool el = false) {
  std::cout << optcstr << "[";
  bool notfirst = false;
  for (auto elem : coll) {
    if (notfirst) {
      std::cout << ',' << elem;
    } else {
      std::cout << elem;
      notfirst = true;
    }
  }
  std::cout << "]";
  if (el)
    std::cout << std::endl;
}


namespace psi{ namespace so_mp2{

extern "C" PSI_API
int read_options(std::string name, Options &options)
{
    if (name == "SO_MP2" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}


extern "C" PSI_API
SharedWavefunction so_mp2(SharedWavefunction ref_wfn, Options& options)
{
    /*
     * This plugin shows a simple way of obtaining MO basis integrals, directly from a DPD buffer.  It is also
     * possible to generate integrals with labels (IWL) formatted files, but that's not shown here.
     */
    int print = options.get_int("PRINT");

    // Grab the global (default) PSIO object, for file I/O
    std::shared_ptr<PSIO> psio(_default_psio_lib_);

    // Have the reference (SCF) wavefunction, ref_wfn
    if(!ref_wfn) throw PSIEXCEPTION("SCF has not been run yet!");

    // Quickly check that there are no open shell orbitals here...
    int nirrep  = ref_wfn->nirrep();

    // For now, we'll just transform for closed shells and generate all integrals.  For more elaborate use of the
    // LibTrans object, check out the plugin_mp2 example in the test suite.
    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(ref_wfn, spaces, IntegralTransform::TransformationType::Restricted);
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());

    /*
     * Now, loop over the DPD buffer, printing the integrals
     */
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"),
                  ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");

    size_t nmo = ref_wfn->nmo();
    size_t nmo4 = nmo * nmo *nmo *nmo;
    std::vector<double> mo_ints(nmo4, 0.0);

    //fancy inline defining of wavefunction! :D
    auto four_idx = [&](size_t p, size_t q, size_t r, size_t s, size_t dim) -> size_t
    {
      size_t dim2 = dim * dim;
      size_t dim3 = dim2 * dim;

      return (p*dim3 + q*dim2 + r*dim + s);
    };

    for(int h = 0; h < nirrep; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            int p = K.params->roworb[h][pq][0];
            int q = K.params->roworb[h][pq][1];
            int psym = K.params->psym[p];
            int qsym = K.params->qsym[q];
            int prel = p - K.params->poff[psym];
            int qrel = q - K.params->qoff[qsym];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                int r = K.params->colorb[h][rs][0];
                int s = K.params->colorb[h][rs][1];
                int rsym = K.params->rsym[r];
                int ssym = K.params->ssym[s];
                int rrel = r - K.params->roff[rsym];
                int srel = s - K.params->soff[ssym];
                // Print out the absolute orbital numbers, the relative (within irrep)
                // numbers, the symmetries, and the integral itself
                // psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f, "
                //                  "symmetries = (%1d %1d | %1d %1d), "
                //                  "relative indices = (%2d %2d | %2d %2d)\n",
                //                  p, q, r, s, K.matrix[h][pq][rs],
                //                  psym, qsym, rsym, ssym,
                //                  prel, qrel, rrel, srel);



                mo_ints[four_idx(p,q,r,s,nmo)] = K.matrix[h][pq][rs];

            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    PRINT_ELEMENTS(mo_ints);

    size_t nso = 2*nmo;
    std::vector<std::pair<size_t, int> > so_lables(nso);
    for(size_t n=0; n<nmo; n++){
      so_lables[2 * n] = std::make_pair(n, 0);
      so_lables[2*n + 1] = std::make_pair(n, 1);
    }
    //now make actual integrals vector!
    size_t nso4 = nso*nso*nso*nso;
    std::vector<double> so_ints(nso4, 0.0);
    for(size_t p=0; p < nso; p++){
      size_t p_orb = so_lables[p].first;
      int p_sp = so_lables[p].second;
      for(size_t q=0; q < nso; q++){
        size_t q_orb = so_lables[q].first;
        int q_sp = so_lables[q].second;
        for(size_t r=0; r < nso; r++){
          size_t r_orb = so_lables[r].first;
          int r_sp = so_lables[r].second;
          for(size_t s=0; s < nso; s++){
            size_t s_orb = so_lables[s].first;
            int s_sp = so_lables[s].second;
            double integral = 0.0;

            if( (p_sp==r_sp) and (q_sp==s_sp) ){
              integral += mo_ints[four_idx(p_orb, r_orb, q_orb, s_orb, nmo)];
            }
            if( (p_sp==s_sp) and (q_sp==r_sp) ){
              integral -= mo_ints[four_idx(p_orb, s_orb, q_orb, r_orb, nmo)];
            }
            so_ints[four_idx(p,q,r,s,nso)] = integral;
            //std::cout << "numerator : " << num << std::endl;
          }
        }
      }
    }

    //PRINT_ELEMENTS(so_ints);

    std::vector<double> eps(nso);
    SharedVector eps_mo = ref_wfn->epsilon_a();

    for(size_t i = 0; i < nso; i++){
      size_t i_orb = so_lables[i].first;
      eps[i] = eps_mo->get(i_orb);
    }

    PRINT_ELEMENTS(eps);

    int na = ref_wfn->nalpha();
    int nb = ref_wfn->nbeta();
    int nocc = na + nb;
    // Assume RHF singet

    std::vector<int> occ;
    std::vector<int> vir;

    INSERT_ELEMENTS(occ, 0, nocc);
    INSERT_ELEMENTS(vir, nocc, nso);

    PRINT_ELEMENTS(occ);
    PRINT_ELEMENTS(vir);


    double Emp2 = 0.0;

    for(int int_i: occ){
      for(int int_j: occ){
        for(int int_a: vir){
          for(int int_b: vir){
            double num = so_ints[four_idx(int_i ,int_j ,int_a, int_b, nso)];
            //std::cout << "numerator : " << num << std::endl;
            double denom = (eps[int_i] + eps[int_j] - eps[int_a] - eps[int_b]);

            Emp2 += 0.25 * num * num /denom;
          }
        }
      }
    }

    double rhf = ref_wfn->reference_energy();
    outfile->Printf("\n    RHF energy     %1.14f",rhf);
    outfile->Printf("\n    MP2 energy     %1.14f",Emp2);
    outfile->Printf("\n    BOTH           %1.14f",rhf + Emp2);

    Process::environment.globals["CURRENT ENERGY"] = rhf + Emp2;

    return ref_wfn;
}

}} // End Namespaces
