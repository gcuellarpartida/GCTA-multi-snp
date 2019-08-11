#include "gcta.h"

void gcta::read_metafile(string metafile, bool GC, double GC_val) {
    cout << "\nReading GWAS summary-level statistics from [" + metafile + "] ..." << endl;
    ifstream Meta(metafile.c_str());
    if (!Meta) throw ("Error: can not open the file [" + metafile + "] to read.");

    int i = 0, count = 0;
    double f_buf = 0.0, b_buf = 0.0, se_buf = 0.0, p_buf = 0.0, N_buf = 0.0, Vp_buf = 0.0, GC_buf = 0.0, chi_buf = 0.0, h_buf = 0.0;
    string A1_buf, A2_buf;
    string snp_buf, str_buf0, str_buf;

    vector<string> snplist, vs_buf, bad_snp;
    vector<string> ref_A_buf, bad_A1, bad_A2, bad_refA;
    vector<double> freq_buf, beta_buf, beta_se_buf, pval_buf, N_o_buf, Vp_v_buf, GC_v_buf;
    map<string, int>::iterator iter;
    getline(Meta, str_buf); // the header line
    if (StrFunc::split_string(str_buf, vs_buf) < 7) throw ("Error: format error in the input file [" + metafile + "].");
    _jma_Vp = 0.0;
    _GC_val = -1;
    while (Meta) {
        getline(Meta, str_buf0);
        stringstream iss(str_buf0);
        iss >> snp_buf >> A1_buf >> A2_buf;
        StrFunc::to_upper(A1_buf);
        StrFunc::to_upper(A2_buf);
        iss >> str_buf;
        f_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == ".") continue;
        b_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == "." || str_buf == "0") continue;
        se_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == ".") continue;
        p_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == ".") continue;
        N_buf = atof(str_buf.c_str());
        if (N_buf < 10) throw ("Error: invalid sample size in line:\n\"" + str_buf0 + "\"");
        if (Meta.eof()) break;
        iter = _snp_name_map.find(snp_buf);
        h_buf = 2.0 * f_buf * (1.0 - f_buf);
        Vp_buf = h_buf * N_buf * se_buf * se_buf + h_buf * b_buf * b_buf * N_buf / (N_buf - 1.0);
        if (Vp_buf < 0.0) throw ("Error: in line:\n\"" + str_buf0 + "\"");
        Vp_v_buf.push_back(Vp_buf);
        if (GC) {
            GC_buf = b_buf * b_buf / se_buf / se_buf;
            if (GC_buf < 0) throw ("Error: in line:\n\"" + str_buf0 + "\"");
            GC_v_buf.push_back(GC_buf);
        }
        count++;
        if (iter == _snp_name_map.end()) continue;
        i = iter->second;
        if (A1_buf != _allele1[i] && A1_buf != _allele2[i]) {
            bad_snp.push_back(_snp_name[i]);
            bad_A1.push_back(_allele1[i]);
            bad_A2.push_back(_allele2[i]);
            bad_refA.push_back(A1_buf);
            continue;
        }
        snplist.push_back(snp_buf);
        ref_A_buf.push_back(A1_buf);
        freq_buf.push_back(f_buf);
        beta_buf.push_back(b_buf);
        beta_se_buf.push_back(se_buf);
        pval_buf.push_back(p_buf);
        N_o_buf.push_back(N_buf);
    }
    Meta.close();
    cout << "GWAS summary statistics of " << count << " SNPs read from [" + metafile + "]." << endl;
    _jma_Vp = CommFunc::median(Vp_v_buf);
    cout << "Phenotypic variance estimated from summary statistics of all " << count << " SNPs: " << _jma_Vp << " (variance of logit for case-control studies)." << endl;
    if (GC) {
        if (GC_val > 0) {
            _GC_val = GC_val;
            cout << "User specified genomic inflation factor: " << _GC_val << endl;
        } else {
            _GC_val = CommFunc::median(GC_v_buf) / 0.455;
            cout << "Genomic inflation factor calcualted from " << count << " SNPs: " << _GC_val << endl;
        }
        cout << "p-values wil be adjusted by the genomic control approach." << endl;
    }

    cout << "Matching the GWAS meta-analysis results to the genotype data ..." << endl;
    update_id_map_kp(snplist, _snp_name_map, _include);
    vector<int> indx(_include.size());
    map<string, int> id_map;
    for (i = 0; i < snplist.size(); i++) id_map.insert(pair<string, int>(snplist[i], i));
    for (i = 0; i < _include.size(); i++) {
        iter = id_map.find(_snp_name[_include[i]]);
        indx[i] = iter->second;
        _ref_A[_include[i]] = ref_A_buf[iter->second];
        if (!_mu.empty() && ref_A_buf[iter->second] == _allele2[_include[i]]) _mu[_include[i]] = 2.0 - _mu[_include[i]];
    }
    if (!bad_snp.empty()) {
        string badsnpfile = _out + ".badsnps";
        ofstream obadsnp(badsnpfile.c_str());
        obadsnp << "SNP\tA1\tA2\tRefA" << endl;
        for (i = 0; i < bad_snp.size(); i++) obadsnp << bad_snp[i] << "\t" << bad_A1[i] << "\t" << bad_A2[i] << "\t" << bad_refA[i] << endl;
        obadsnp.close();
        cout << "Warning: can't match the reference alleles of " << bad_snp.size() << " SNPs to those in the genotype data. These SNPs have been saved in [" + badsnpfile + "]." << endl;
    }
    if (_include.empty()) throw ("Error: all the SNPs in the GWAS meta-analysis results can't be found in the genotype data.");
    else cout << _include.size() << " SNPs are matched to the genotype data." << endl;

    if (_mu.empty()) calcu_mu();

    _freq.resize(_include.size());
    _beta.resize(_include.size());
    _beta_se.resize(_include.size());
    _pval.resize(_include.size());
    _N_o.resize(_include.size());
    for (i = 0; i < _include.size(); i++) {
        _freq[i] = freq_buf[indx[i]];
        _beta[i] = beta_buf[indx[i]];
        _beta_se[i] = beta_se_buf[indx[i]];
        if (GC) {
            chi_buf = _beta[i] / _beta_se[i];
            _pval[i] = StatFunc::pchisq(chi_buf * chi_buf / _GC_val, 1);
        } else _pval[i] = pval_buf[indx[i]];
        _N_o[i] = N_o_buf[indx[i]];
    }
    _jma_Ve = _jma_Vp;
}

void gcta::init_massoc(string metafile, bool GC, double GC_val)
{
    read_metafile(metafile, GC, GC_val);

    int i = 0, j = 0, n = _keep.size(), m = _include.size();
    cout << "Calculating the variance of SNP genotypes ..." << endl;
    _MSX_B.resize(m);
    _Nd.resize(m);

    if (_mu.empty()) calcu_mu();
    #pragma omp parallel for
    for (i = 0; i < m; i++){
        eigenVector x;
        makex_eigenVector(i, x, true, true);
        _MSX_B[i] = x.squaredNorm() / (double)n;
    }
    if (_jma_actual_geno) {
        _MSX = _MSX_B;
        _Nd = _N_o;
    } else {
        _MSX = 2.0 * _freq.array()*(1.0 - _freq.array());
        for (i = 0; i < m; i++) _Nd[i] = (_jma_Vp - _MSX[i] * _beta[i] * _beta[i]) / (_MSX[i] * _beta_se[i] * _beta_se[i]) + 1; // revised by JY 25/11/13 according to Eq. 13, Yang et al. 2012 NG
    }
}

void gcta::read_fixed_snp(string snplistfile, string msg, vector<int> &pgiven, vector<int> &remain) {
    int i = 0, j = 0;
    vector<string> givenSNPs;
    read_snplist(snplistfile, givenSNPs, msg);
    if (givenSNPs.empty()) throw ("Error: failed to read any SNP from the file [" + snplistfile + "].");
    map<string, int> givenSNPs_map;
    pgiven.clear();
    remain.clear();
    for (i = 0; i < givenSNPs.size(); i++) givenSNPs_map.insert(pair<string, int>(givenSNPs[i], i));
    for (i = 0; i < _include.size(); i++) {
        if (givenSNPs_map.find(_snp_name[_include[i]]) != givenSNPs_map.end()) pgiven.push_back(i);
        else remain.push_back(i);
    }
    if (pgiven.size() > 0) cout << pgiven.size() << " of them are matched to the genotype and summary data." << endl;
    else throw ("Error: none of the given SNPs can be matched to the genotype and summary data.");
}

void gcta::run_massoc_slct(string metafile, int wind_size, double p_cutoff, double collinear, int top_SNPs, bool joint_only, bool GC, double GC_val, bool actual_geno, int mld_slct_alg)
{
    _jma_actual_geno = actual_geno;
    _jma_wind_size = wind_size;
    _jma_p_cutoff = p_cutoff;
    _jma_collinear = collinear;
    _jma_snpnum_backward = 0;
    _jma_snpnum_collienar = 0;
    init_massoc(metafile, GC, GC_val);
    if (top_SNPs < 0) top_SNPs = 1e30;
    else {
        _jma_p_cutoff = 0.5;
        cout << "The threshold p-value has been set to 0.5 because of the --massoc-top-SNPs option." << endl;
    }
    int i = 0, j = 0;
    vector<int> slct, remain;
    eigenVector bC, bC_se, pC;
    cout << endl;
    if (!joint_only && mld_slct_alg<2) {
        cout << "Performing "<< ((mld_slct_alg==0)?"stepwise":"forward") << " model selection on " << _include.size() << " SNPs to select association signals ... (p cutoff = " << _jma_p_cutoff << "; ";
        cout << "collinearity cutoff = " << _jma_collinear << ")"<< endl;
        if (!_jma_actual_geno) cout << "(Assuming complete linkage equilibrium between SNPs which are more than " << _jma_wind_size / 1e6 << "Mb away from each other)" << endl;
        stepwise_slct(slct, remain, bC, bC_se, pC, mld_slct_alg, top_SNPs);
        if (slct.empty()) {
            cout << "No SNPs have been selected." << endl;
            return;
        }
    } else { // mld_slct_alg = 2 for backward selection
        for (i = 0; i < _include.size(); i++) slct.push_back(i);
        if (mld_slct_alg==2) {
            cout << "Performing backward selection on " << _include.size() << " SNPs at threshold p-value = " << _jma_p_cutoff << " ..." << endl;
            slct_stay(slct, bC, bC_se, pC);
        }
    }

    // joint analysis
    eigenVector bJ, bJ_se, pJ;
    cout << "Performing joint analysis on all the " << slct.size();
    if (joint_only) cout << " SNPs ..." << endl;
    else cout << " selected signals ..." << endl;
    if (slct.size() >= _keep.size()) throw ("Error: too many SNPs. The number of SNPs in a joint analysis should not be larger than the sample size.");
    massoc_joint(slct, bJ, bJ_se, pJ);
    eigenMatrix rval(slct.size(), slct.size());
    LD_rval(slct, rval);
    if (_jma_actual_geno) cout << "Residual varaince = " << _jma_Ve << endl;
    massoc_slct_output(joint_only, slct, bJ, bJ_se, pJ, rval);

    // output conditional results
    if (!joint_only && !mld_slct_alg==2) {
        massoc_cond_output(remain, bC, bC_se, pC);
        cout << "(" << _jma_snpnum_backward << " SNPs eliminated by backward selection and " << _jma_snpnum_collienar << " SNPs filtered by collinearity test are not included in the output)" << endl;
    }
}

void gcta::run_massoc_cond(string metafile, string snplistfile, int wind_size, double collinear, bool GC, double GC_val, bool acutal_geno) {
    _jma_actual_geno = acutal_geno;
    _jma_wind_size = wind_size;
    _jma_collinear = collinear;
    init_massoc(metafile, GC, GC_val);
    vector<int> pgiven, remain;
    read_fixed_snp(snplistfile, "given SNPs", pgiven, remain);

    eigenVector bC, bC_se, pC;
    cout << "Performing single-SNP association analysis conditional on the " << pgiven.size() << " given SNPs ... ";
    cout << "(collinearity cutoff = " << _jma_collinear << ")" << endl;
    if (!_jma_actual_geno) cout << "(Assuming complete linkage equilibrium between SNPs which are more than " << _jma_wind_size / 1e6 << "Mb away from each other)" << endl;

    string filename = _out + ".given.cojo";
    cout<<"Saving the summary statistics of the given SNPs in the file [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) throw ("Can not open the file [" + filename + "] to write.");
    ofile << "Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp" << endl;
    int i = 0, j = 0;
    for (i = 0; i < pgiven.size(); i++) {
        j = pgiven[i];
        ofile << _chr[_include[j]] << "\t" << _snp_name[_include[j]] << "\t" << _bp[_include[j]] << "\t" << _ref_A[_include[j]] << "\t" << _freq[j] << "\t" << _beta[j] << "\t" << _beta_se[j] << "\t" << _pval[j] << endl;
    }
    ofile.close();
    
    massoc_cond(pgiven, remain, bC, bC_se, pC);
    massoc_cond_output(remain, bC, bC_se, pC);
}

void gcta::massoc_slct_output(bool joint_only, vector<int> &slct, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, eigenMatrix &rval)
{
    string filename = _out + ".jma.cojo";
    if (joint_only) cout << "Saving the joint analysis result of " << slct.size() << " SNPs to [" + filename + "] ..." << endl;
    else cout << "Saving the " << slct.size() << " independent signals to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) throw ("Can not open the file [" + filename + "] to write.");
    ofile << "Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp" << ((_GC_val > 0) ? "_GC" : "") << "\tn\tfreq_geno\tbJ\tbJ_se\tpJ" << ((_GC_val > 0) ? "_GC" : "") << "\tLD_r" << endl;
    int i = 0, j = 0;
    for (i = 0; i < slct.size(); i++) {
        j = slct[i];
        ofile << _chr[_include[j]] << "\t" << _snp_name[_include[j]] << "\t" << _bp[_include[j]] << "\t";
        ofile << _ref_A[_include[j]] << "\t" << _freq[j] << "\t" << _beta[j] << "\t" << _beta_se[j] << "\t";
        ofile << _pval[j] << "\t" << _Nd[j] << "\t" << 0.5 * _mu[_include[j]] << "\t" << bJ[i] << "\t" << bJ_se[i] << "\t" << pJ[i] << "\t";
        if (i == slct.size() - 1) ofile << 0 << endl;
        else ofile << rval(i, i + 1) << endl;
    }
    ofile.close();

    filename = _out + ".ldr.cojo";
    if (joint_only) cout << "Saving the LD structure of " << slct.size() << " SNPs to [" + filename + "] ..." << endl;
    else cout << "Saving the LD structure of " << slct.size() << " independent signals to [" + filename + "] ..." << endl;
    ofile.open(filename.c_str());
    if (!ofile) throw ("Can not open the file [" + filename + "] to write.");
    ofile << "SNP\t";
    for (i = 0; i < slct.size(); i++) ofile << _snp_name[_include[slct[i]]] << "\t";
    ofile << endl;
    for (i = 0; i < slct.size(); i++) {
        ofile << _snp_name[_include[slct[i]]] << "\t";
        for (j = 0; j < slct.size(); j++) ofile << rval(i, j) << "\t";
        ofile << endl;
    }
    ofile.close();
}

void gcta::massoc_cond_output(vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC)
{
    int i = 0, j = 0;
    string filename = _out + ".cma.cojo";
    cout << "Saving the conditional analysis results of " << remain.size() << " remaining SNPs to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) throw ("Can not open the file [" + filename + "] to write.");
    ofile << "Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp" << ((_GC_val > 0) ? "_GC" : "") << "\tn\tfreq_geno\tbC\tbC_se\tpC" << ((_GC_val > 0) ? "_GC" : "") << endl;
    for (i = 0; i < remain.size(); i++) {
        j = remain[i];
        ofile << _chr[_include[j]] << "\t" << _snp_name[_include[j]] << "\t" << _bp[_include[j]] << "\t";
        ofile << _ref_A[_include[j]] << "\t" << _freq[j] << "\t" << _beta[j] << "\t" << _beta_se[j] << "\t";
        ofile << _pval[j] << "\t" << _Nd[j] << "\t" << 0.5 * _mu[_include[j]] << "\t";
        if (pC[i] > 1.5) ofile << "NA\tNA\tNA" << endl;
        else ofile << bC[i] << "\t" << bC_se[i] << "\t" << pC[i] << endl;
    }
    ofile.close();
}

void gcta::stepwise_slct(vector<int> &slct, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, int mld_slct_alg, int top_SNPs)
{
    int i = 0, i_buf = 0;
    vector<double> p_buf;
    eigenVector2Vector(_pval, p_buf);
    int m = min_element(p_buf.begin(), p_buf.end()) - p_buf.begin();
    if (p_buf[m] >= _jma_p_cutoff) return;
    slct.push_back(m);
    for (i = 0; i < _include.size(); i++) {
        if (i != m) remain.push_back(i);
    }
    int prev_num = 0;
    if (mld_slct_alg==0 && _jma_p_cutoff > 1e-3){
        cout << "Switched to perform a forward model selection because the significance level is too low..." << endl;
        mld_slct_alg=1;
    }
    while (!remain.empty()) {
        if (slct_entry(slct, remain, bC, bC_se, pC)) {
            if (mld_slct_alg==0) slct_stay(slct, bC, bC_se, pC);
        }
        else break;
        if (slct.size() % 5 == 0 && slct.size() > prev_num) cout << slct.size() << " associated SNPs have been selected." << endl;
        if (slct.size() > prev_num) prev_num = slct.size();
        if (slct.size() >= top_SNPs) break;
    }
    if (_jma_p_cutoff > 1e-3) {
        cout << "Performing backward elimination..." << endl;
        slct_stay(slct, bC, bC_se, pC);
    }
    cout << "Finally, " << slct.size() << " associated SNPs are selected." << endl;
}

bool gcta::slct_entry(vector<int> &slct, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC) {
    int i = 0, m = 0;
    massoc_cond(slct, remain, bC, bC_se, pC);
    vector<double> pC_buf;
    eigenVector2Vector(pC, pC_buf);
    while (true) {
        m = min_element(pC_buf.begin(), pC_buf.end()) - pC_buf.begin();
        if (pC_buf[m] >= _jma_p_cutoff) return (false);
        if (insert_B_and_Z(slct, remain[m])) {
            slct.push_back(remain[m]);
            stable_sort(slct.begin(), slct.end());
            remain.erase(remain.begin() + m);
            return (true);
        }
        pC_buf.erase(pC_buf.begin() + m);
        remain.erase(remain.begin() + m);
    }
}

void gcta::slct_stay(vector<int> &slct, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ) {
    if (_B_N.cols() < 1) {
        if (!init_B(slct)) throw ("Error: there is a collinearity problem of the given list of SNPs.\nYou can try the option --massoc-slct to get rid of one of each pair of highly correlated SNPs.");
    }

    vector<double> pJ_buf;
    while (!slct.empty()) {
        massoc_joint(slct, bJ, bJ_se, pJ);
        eigenVector2Vector(pJ, pJ_buf);
        int m = max_element(pJ_buf.begin(), pJ_buf.end()) - pJ_buf.begin();
        if (pJ[m] > _jma_p_cutoff) {
            _jma_snpnum_backward++;
            erase_B_and_Z(slct, slct[m]);
            slct.erase(slct.begin() + m);
        } else break;
    }
}

void gcta::eigenVector2Vector(eigenVector &x, vector<double> &y) {
    y.resize(x.size());
    for (int i = 0; i < x.size(); i++) y[i] = x[i];
}

double gcta::massoc_calcu_Ve(const vector<int> &slct, eigenVector &bJ, eigenVector &b) {
    double Ve = 0.0;
    int n = bJ.size();
    vector<double> Nd_buf(n);
    for (int k = 0; k < n; k++) {
        Nd_buf[k] = _Nd[slct[k]];
        Ve += _D_N[k] * bJ[k] * b[k];
    }
    double d_buf = CommFunc::median(Nd_buf);
    if (d_buf - n < 1) throw ("Error: no degree of freedom left for the residues, the model is over-fitting. Please specify a more stringent p cutoff value.");
    Ve = ((d_buf - 1) * _jma_Vp - Ve) / (d_buf - n);
    if (Ve <= 0.0) throw ("Error: residual variance is out of boundary, the model is over-fitting. Please specify a more stringent p cutoff value.");
    return Ve;
}

void gcta::massoc_joint(const vector<int> &indx, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ) {
    if (_B_N.cols() < 1) {
        if (!init_B(indx)) throw ("Error: there is a collinearity problem of the given list of SNPs.\nYou can try the option --massoc-slct to get rid of one of each pair of highly correlated SNPs.");
    }

    int i = 0, n = indx.size();
    double chisq = 0.0;
    eigenVector b(n);
    for (i = 0; i < n; i++) b[i] = _beta[indx[i]];
    bJ.resize(n);
    bJ_se.resize(n);
    pJ.resize(n);
    bJ = _B_N_i * _D_N.asDiagonal() * b;
    bJ_se = _B_N_i.diagonal();
    pJ = eigenVector::Ones(n);
    if (_jma_actual_geno) _jma_Ve = massoc_calcu_Ve(indx, bJ, b);
    bJ_se *= _jma_Ve;
    for (i = 0; i < n; i++) {
        if (bJ_se[i] > 1.0e-7) {
            bJ_se[i] = sqrt(bJ_se[i]);
            chisq = bJ[i] / bJ_se[i];
            if (_GC_val > 0) pJ[i] = StatFunc::pchisq(chisq * chisq / _GC_val, 1);
            else pJ[i] = StatFunc::pchisq(chisq*chisq, 1);
        } else {
            bJ[i] = 0.0;
            bJ_se[i] = 0.0;
        }
    }
}

void gcta::massoc_cond(const vector<int> &slct, const vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC) {
    if (_B_N.cols() < 1) {
        if (!init_B(slct)) throw ("Error: there is a collinearity problem of the given list of SNPs.\nYou can try the option --massoc-slct to get rid of one of each pair of highly correlated SNPs.");
    }
    if (_Z_N.cols() < 1) init_Z(slct);

    int i = 0, j = 0, n = slct.size();
    double chisq = 0.0;
    eigenVector b(n), se(n);
    for (i = 0; i < n; i++) {
        b[i] = _beta[slct[i]];
        se[i] = _beta_se[slct[i]];
    }
    eigenVector bJ1 = _B_N_i * _D_N.asDiagonal() * b;
    if (_jma_actual_geno) _jma_Ve = massoc_calcu_Ve(slct, bJ1, b);

    double B2 = 0.0;
    bC = eigenVector::Zero(remain.size());
    bC_se = eigenVector::Zero(remain.size());
    pC = eigenVector::Constant(remain.size(), 2);
    eigenVector Z_Bi(n), Z_Bi_buf(n);
    for (i = 0; i < remain.size(); i++) {
        j = remain[i];
        B2 = _MSX[j] * _Nd[j];
        if (!CommFunc::FloatEqual(B2, 0.0)) {
            Z_Bi = _Z_N.col(j).transpose() * _B_N_i;
            Z_Bi_buf = _Z.col(j).transpose() * _B_i;
            if (_Z.col(j).dot(Z_Bi_buf) / _MSX_B[j] < _jma_collinear) {
                bC[i] = _beta[j] - Z_Bi.cwiseProduct(_D_N).dot(b) / B2;
                bC_se[i] = (B2 - _Z_N.col(j).dot(Z_Bi)) / (B2 * B2);
            }
        }
        if (_jma_actual_geno) bC_se[i] *= _jma_Ve - (B2 * bC[i] * _beta[j]) / (_Nd[j] - n - 1);
        else bC_se[i] *= _jma_Ve;
        if (bC_se[i] > 1e-7) {
            bC_se[i] = sqrt(bC_se[i]);
            chisq = bC[i] / bC_se[i];
            if (_GC_val > 0) pC[i] = StatFunc::pchisq(chisq * chisq / _GC_val, 1);
            else pC[i] = StatFunc::pchisq(chisq*chisq, 1);
        }
    }
}

bool gcta::init_B(const vector<int> &indx)
{
    int i = 0, j = 0, k = 0, n = _keep.size();
    double d_buf = 0.0;
    eigenVector diagB(indx.size());
    _B.resize(indx.size(), indx.size());
    _B_N.resize(indx.size(), indx.size());
    _D_N.resize(indx.size());
    eigenVector x_i(_keep.size()), x_j(_keep.size());
    for (i = 0; i < indx.size(); i++) {
        _D_N[i] = _MSX[indx[i]] * _Nd[indx[i]];
        _B.startVec(i);
        _B_N.startVec(i);
        _B.insertBack(i, i) = _MSX_B[indx[i]];
        _B_N.insertBack(i, i) = _D_N[i];
        diagB[i] = _MSX_B[indx[i]];
        makex_eigenVector(indx[i], x_i, false, true);
        for (j = i + 1; j < indx.size(); j++) {
            if (_jma_actual_geno || (_chr[_include[indx[i]]] == _chr[_include[indx[j]]] && abs(_bp[_include[indx[i]]] - _bp[_include[indx[j]]]) < _jma_wind_size)) {
                makex_eigenVector(indx[j], x_j, false, true);
                d_buf = x_i.dot(x_j) / (double)n;
                _B.insertBack(j, i) = d_buf;
                _B_N.insertBack(j, i) = d_buf * min(_Nd[indx[i]], _Nd[indx[j]]) * sqrt(_MSX[indx[i]] * _MSX[indx[j]] / (_MSX_B[indx[i]] * _MSX_B[indx[j]]));
            }
        }
    }
    _B.finalize();
    _B_N.finalize();

    LDLT<eigenMatrix> ldlt_B(_B);
    if (ldlt_B.vectorD().minCoeff() < 0 || sqrt(ldlt_B.vectorD().maxCoeff() / ldlt_B.vectorD().minCoeff()) > 30) return false;
    _B_i = eigenMatrix::Identity(indx.size(), indx.size());
    ldlt_B.solveInPlace(_B_i);
    if ((1 - eigenVector::Constant(indx.size(), 1).array() / (diagB.array() * _B_i.diagonal().array())).maxCoeff() > _jma_collinear) return false;
    LDLT<eigenMatrix> ldlt_B_N(_B_N);
    _B_N_i = eigenMatrix::Identity(indx.size(), indx.size());
    ldlt_B_N.solveInPlace(_B_N_i);
    return true;
}

void gcta::init_Z(const vector<int> &indx)
{
    int i = 0, j = 0, n = _keep.size();
    double d_buf = 0.0;
    _Z.resize(indx.size(), _include.size());
    _Z_N.resize(indx.size(), _include.size());
    eigenVector x_i(_keep.size()), x_j(_keep.size());
    for (j = 0; j < _include.size(); j++) {
        _Z.startVec(j);
        _Z_N.startVec(j);
        makex_eigenVector(j, x_j, false, true);
        for (i = 0; i < indx.size(); i++) {
            if (_jma_actual_geno || (indx[i] != j && _chr[_include[indx[i]]] == _chr[_include[j]] && abs(_bp[_include[indx[i]]] - _bp[_include[j]]) < _jma_wind_size)) {
                makex_eigenVector(indx[i], x_i, false, true);
                d_buf = x_j.dot(x_i) / (double)n;
                _Z.insertBack(i, j) = d_buf;
                _Z_N.insertBack(i, j) = d_buf * min(_Nd[indx[i]], _Nd[j]) * sqrt(_MSX[indx[i]] * _MSX[j] / (_MSX_B[indx[i]] * _MSX_B[j])); // added by Jian Yang 18/12/2013
            }
        }
    }
    _Z.finalize();
    _Z_N.finalize();
}

bool gcta::insert_B_and_Z(const vector<int> &indx, int insert_indx)
{
    int i = 0, j = 0, n = _keep.size();
    double d_buf = 0.0;
    vector<int> ix(indx);
    ix.push_back(insert_indx);
    stable_sort(ix.begin(), ix.end());
    eigenSparseMat B_buf(_B), B_N_buf(_B_N);
    _B.resize(ix.size(), ix.size());
    _B_N.resize(ix.size(), ix.size());
    bool get_insert_col = false, get_insert_row = false;
    int pos = find(ix.begin(), ix.end(), insert_indx) - ix.begin();
    eigenVector diagB(ix.size());
    eigenVector x_i(_keep.size()), x_j(_keep.size());
    for (j = 0; j < ix.size(); j++) {
        _B.startVec(j);
        _B_N.startVec(j);
        _B.insertBack(j, j) = _MSX_B[ix[j]];
        _B_N.insertBack(j, j) = _MSX[ix[j]] * _Nd[ix[j]];
        diagB[j] = _MSX_B[ix[j]];
        if (insert_indx == ix[j]) get_insert_col = true;
        get_insert_row = get_insert_col;
        makex_eigenVector(ix[j], x_j, false, true);
        for (i = j + 1; i < ix.size(); i++) {
            if (insert_indx == ix[i]) get_insert_row = true;
            if (insert_indx == ix[j] || insert_indx == ix[i]) {
                if (_jma_actual_geno || (_chr[_include[ix[i]]] == _chr[_include[ix[j]]] && abs(_bp[_include[ix[i]]] - _bp[_include[ix[j]]]) < _jma_wind_size)) {
                    makex_eigenVector(ix[i], x_i, false, true);
                    d_buf = x_i.dot(x_j) / (double)n;
                    _B.insertBack(i, j) = d_buf;
                    _B_N.insertBack(i, j) = d_buf * min(_Nd[ix[i]], _Nd[ix[j]]) * sqrt(_MSX[ix[i]] * _MSX[ix[j]] / (_MSX_B[ix[i]] * _MSX_B[ix[j]]));
                }
            } else {
                if (B_buf.coeff(i - get_insert_row, j - get_insert_col) != 0) {
                    _B.insertBack(i, j) = B_buf.coeff(i - get_insert_row, j - get_insert_col);
                    _B_N.insertBack(i, j) = B_N_buf.coeff(i - get_insert_row, j - get_insert_col);
                }
            }
        }
    }
    _B.finalize();
    _B_N.finalize();
    LDLT<eigenMatrix> ldlt_B(_B);
    _B_i = eigenMatrix::Identity(ix.size(), ix.size());
    ldlt_B.solveInPlace(_B_i);
    if (ldlt_B.vectorD().minCoeff() < 0 || sqrt(ldlt_B.vectorD().maxCoeff() / ldlt_B.vectorD().minCoeff()) > 30 || (1 - eigenVector::Constant(ix.size(), 1).array() / (diagB.array() * _B_i.diagonal().array())).maxCoeff() > _jma_collinear) {
        _jma_snpnum_collienar++;
        _B = B_buf;
        _B_N = B_N_buf;
        return false;
    }
    LDLT<eigenMatrix> ldlt_B_N(_B_N);
    _B_N_i = eigenMatrix::Identity(ix.size(), ix.size());
    ldlt_B_N.solveInPlace(_B_N_i);
    _D_N.resize(ix.size());
    for (j = 0; j < ix.size(); j++) {
        _D_N[j] = _MSX[ix[j]] * _Nd[ix[j]];
    }

    if (_Z_N.cols() < 1) return true;
    eigenSparseMat Z_buf(_Z), Z_N_buf(_Z_N);
    _Z.resize(ix.size(), _include.size());
    _Z_N.resize(ix.size(), _include.size());
    for (j = 0; j < _include.size(); j++) {
        _Z.startVec(j);
        _Z_N.startVec(j);
        get_insert_row = false;
        makex_eigenVector(j, x_j, false, true);
        for (i = 0; i < ix.size(); i++) {
            if (insert_indx == ix[i]) {
                if (_jma_actual_geno || (ix[i] != j && _chr[_include[ix[i]]] == _chr[_include[j]] && abs(_bp[_include[ix[i]]] - _bp[_include[j]]) < _jma_wind_size)) {
                    makex_eigenVector(ix[i], x_i, false, true);
                    d_buf = x_j.dot(x_i) / (double)n;
                    _Z.insertBack(i, j) = d_buf;
                    _Z_N.insertBack(i, j) = d_buf * min(_Nd[ix[i]], _Nd[j]) * sqrt(_MSX[ix[i]] * _MSX[j] / (_MSX_B[ix[i]] * _MSX_B[j])); // added by Jian Yang 18/12/2013
                }
                get_insert_row = true;
            } else {
                if (Z_buf.coeff(i - get_insert_row, j) != 0) {
                    _Z.insertBack(i, j) = Z_buf.coeff(i - get_insert_row, j);
                    _Z_N.insertBack(i, j) = Z_N_buf.coeff(i - get_insert_row, j);
                }
            }
        }
    }
    _Z.finalize();
    _Z_N.finalize();

    return true;
}

void gcta::erase_B_and_Z(const vector<int> &indx, int erase_indx) {
    int i = 0, j = 0;
    eigenSparseMat B_dense(_B), B_N_dense(_B_N);
    _B.resize(indx.size() - 1, indx.size() - 1);
    _B_N.resize(indx.size() - 1, indx.size() - 1);
    _D_N.resize(indx.size() - 1);
    int pos = find(indx.begin(), indx.end(), erase_indx) - indx.begin();
    bool get_insert_col = false, get_insert_row = false;
    for (j = 0; j < indx.size(); j++) {
        if (erase_indx == indx[j]) {
            get_insert_col = true;
            continue;
        }
        _B.startVec(j - get_insert_col);
        _B_N.startVec(j - get_insert_col);
        _D_N[j - get_insert_col] = _MSX[indx[j]] * _Nd[indx[j]];
        get_insert_row = get_insert_col;
        for (i = j; i < indx.size(); i++) {
            if (erase_indx == indx[i]) {
                get_insert_row = true;
                continue;
            }
            if (B_dense.coeff(i, j) != 0) {
                _B.insertBack(i - get_insert_row, j - get_insert_col) = B_dense.coeff(i, j);
                _B_N.insertBack(i - get_insert_row, j - get_insert_col) = B_N_dense.coeff(i, j);
            }
        }
    }
    _B.finalize();
    _B_N.finalize();

    if (_Z_N.cols() < 1) return;

    LDLT<eigenMatrix> ldlt_B(_B);
    _B_i = eigenMatrix::Identity(indx.size() - 1, indx.size() - 1);
    ldlt_B.solveInPlace(_B_i);
    LDLT<eigenMatrix> ldlt_B_N(_B_N);
    _B_N_i = eigenMatrix::Identity(indx.size() - 1, indx.size() - 1);
    ldlt_B_N.solveInPlace(_B_N_i);

    eigenSparseMat Z_buf(_Z), Z_N_buf(_Z_N);
    _Z.resize(indx.size() - 1, _include.size());
    _Z_N.resize(indx.size() - 1, _include.size());
    for (j = 0; j < _include.size(); j++) {
        _Z.startVec(j);
        _Z_N.startVec(j);
        get_insert_row = false;
        for (i = 0; i < indx.size(); i++) {
            if (erase_indx == indx[i]) {
                get_insert_row = true;
                continue;
            }
            if (Z_buf.coeff(i, j) != 0) {
                _Z.insertBack(i - get_insert_row, j) = Z_buf.coeff(i, j);
                _Z_N.insertBack(i - get_insert_row, j) = Z_N_buf.coeff(i, j);
            }
        }
    }
    _Z.finalize();
    _Z_N.finalize();
}

/*
double gcta::crossprod(vector<float> &x_i, vector<float> &x_j) {
    double prod = 0.0;
    for (int i = 0; i < x_i.size(); i++) prod += x_i[i] * x_j[i];
    return (prod / x_i.size());
}
*/

void gcta::LD_rval(const vector<int> &indx, eigenMatrix &rval) {
    int i = 0, j = 0;
    eigenVector sd(indx.size());
    for (i = 0; i < indx.size(); i++) sd[i] = sqrt(_MSX_B[indx[i]]);
    for (j = 0; j < indx.size(); j++) {
        rval(j, j) = 1.0;
        for (i = j + 1; i < indx.size(); i++) rval(i, j) = rval(j, i) = _B.coeff(i, j) / sd[i] / sd[j];
    }
}

void gcta::run_massoc_sblup(string metafile, int wind_size, double lambda)
{
    _jma_wind_size = wind_size;
    init_massoc(metafile, false, -1);

    int j = 0;
    eigenVector bJ;
    cout << "\nPerforming joint analysis on all the " << _include.size() << " SNPs ..." << endl;
    cout << "(Assuming complete linkage equilibrium between SNPs which are more than " << _jma_wind_size / 1e6 << "Mb away from each other)" << endl;
    if (massoc_sblup(lambda, bJ)) {
        string filename = _out + ".sblup.cojo";
        cout << "Saving the joint analysis result of " << _include.size() << " SNPs to [" + filename + "] ..." << endl;
        ofstream ofile(filename.c_str());
        if (!ofile) throw ("Can not open the file [" + filename + "] to write.");
        for (j = 0; j < _include.size(); j++) {
            ofile << _snp_name[_include[j]] << "\t" << _ref_A[_include[j]] << "\t" << _beta[j] << "\t" << bJ[j] << endl;
        }
        ofile.close();
    } else throw ("Jacobi iteration not converged. You can increase the maximum number of iterations by the option --massoc-sblup-maxit.");
}

bool gcta::massoc_sblup(double lambda, eigenVector &bJ)
{
    int i = 0, j = 0, k = 0, n = _keep.size(), m = _include.size();
    double prod = 0.0;
    eigenVector D(_include.size());
    eigenSparseMat B(_include.size(), _include.size());

    vector<int> nz_i(m), nz_j(m);
    for (i = 0; i < _include.size(); i++){
        nz_i[i]=0;
        nz_j[i]=0;
    }
    cout << "Calculating the LD correlation matrix of all the " << _include.size() << " SNPs..." << endl;
    for (i = 0; i < _include.size(); i++) {
        D[i] = _MSX[i] * _Nd[i];
        B.startVec(i);
        B.insertBack(i, i) = D[i] + lambda;
        for (j = i + 1; j < _include.size(); j++) {
            if (_chr[_include[i]] == _chr[_include[j]] && abs(_bp[_include[i]] - _bp[_include[j]]) < _jma_wind_size) {
                B.insertBack(j, i) = 1.0;
            }
        }
    }
    B.finalize();

    eigenVector x_i(_keep.size()), x_j(_keep.size());
    for (i = 0; i < _include.size(); i++) {
        makex_eigenVector(i, x_i, false, true);
        for (j = i + 1; j < _include.size(); j++) {
            if (_chr[_include[i]] == _chr[_include[j]] && abs(_bp[_include[i]] - _bp[_include[j]]) < _jma_wind_size) {
                makex_eigenVector(j, x_j, false, true);
                prod = x_i.dot(x_j) / (double)n;
                B.coeffRef(j, i) = prod * min(_Nd[i], _Nd[j]) * sqrt(_MSX[i] * _MSX[j] / (_MSX_B[i] * _MSX_B[j]));
            }
        }
        if((i + 1) % 1000 == 0 || (i + 1) == _include.size()) cout << i + 1 << " of " << _include.size() << " SNPs.\r";
    }

    cout << "Estimating the joint effects of all SNPs ..." << endl;
    eigenVector Xty = D.array() * _beta.array();
    SimplicialLDLT<eigenSparseMat> solver;
    solver.compute(B);
    if(solver.info()!=Success) {
        throw("Error: decomposition failed! The SNP correlation matrix is not positive definite.");
    }
    bJ = solver.solve(Xty);
    if(solver.info()!=Success) {
        throw("Error: solving failed! Unable to solve the BLUP equation.");
    }

    // debug
    //cout<<"_MSX: "<<_MSX<<endl;
    //cout<<"Nd: "<<_Nd<<endl;
    //cout<<"B: "<<B<<endl;
    //cout<<"D: "<<D<<endl;
    //cout<<"_beta: "<<_beta<<endl;
    //cout<<"bJ: "<<bJ<<endl;

    //cout<<"B inverse: "<<ldlt.solve(eigenMatrix::Identity(_include.size(), _include.size()))<<endl;


    return (true);
}



