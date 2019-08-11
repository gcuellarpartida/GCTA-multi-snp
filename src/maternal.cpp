#include <unordered_map>

#include "gcta.h"

typedef unordered_map<string, vector<int>> parent_map_t;



std::vector<float> gcta::calc_parental_prob(enum SNPGenotype p1_geno, enum SNPGenotype p2_geno, float maf)
{
    std::vector<float> result;
    float p = maf;
    float q = 1 - maf;

    if (p1_geno == INVALID || p2_geno == INVALID) {
        return result;
    }
    if (p1_geno == AA && p2_geno == AA) 
    {
        result = {(2 * p)/(p + 1), (1 - p)/(1 + p), 0};
    }  
    else if ((p1_geno == AA &&  p2_geno == Aa) || (p1_geno == Aa && p2_geno == AA))
    {
        result = {p/(p+1), 1/(p+1), 0}; 
    }
    else if ((p1_geno == AA && p2_geno == aa) || (p1_geno == aa && p2_geno == AA))
    {
        result = {0, 1, 0};   
    } else if (p1_geno == Aa && p2_geno == Aa)
    {
        float denom = 1+p-p*p;
        result = {(p * (1 - 0.5f*p))/denom, 0.5f/denom, 0.5f*(1 - p*p)/denom};
    } else if ((p1_geno == Aa && p2_geno == aa) || (p1_geno == aa && p2_geno == Aa))
    {
        result = {0, 1/(q + 1), q/(q+1)};
    } else if ((p1_geno == aa && p2_geno == aa)) 
    {
        result = {0, (1-q)/(q+1), 2*q/(q+1)};
    } else 
    {
        throw("Unknown Genotype, should never reach here");
    }
    return result;
}

void gcta::calc_parental_prob_snps(int p1, int p2, eigenMatrix &probs)
{
    probs.resize(1, _snp_num);
    int p1_geno = 0, p2_geno = 0;

    int snp_i;
    #pragma omp parallel for private(snp_i)
    for (snp_i = 0; snp_i < _snp_num; snp_i++) 
    {
        float maf = _maf[snp_i];
        p1_geno = 0, p2_geno = 0;
        // 11 = AA, 00 == aa, 01 == Aa, 10 == INVALID  
        if (!_snp_1[snp_i][p1]) {
            p1_geno |= (1 << 0);        
        } 
        if (!_snp_2[snp_i][p1]) {
            p1_geno |= (1 << 1);
        }

        if (!_snp_1[snp_i][p2]) {
            p2_geno |= (1 << 0);        
        } 
        if (!_snp_2[snp_i][p2]) {
            p2_geno |= (1 << 1);
        }
        // calculate probabilities
        std::vector<float> result = calc_parental_prob((enum SNPGenotype) p1_geno, 
            (enum SNPGenotype) p2_geno, 1-maf);
        probs(0, snp_i) = result.size() ? result[1] + 2 * result[2] : 0.0f;
    }    
}

void gcta::dump_genotypes(string filename, string snp)
{
    ofstream outfile;
    outfile.open(filename);
    if (_snp_name_map.find(snp) == _snp_name_map.end()) {
        throw ("No such SNP found!");
    } 
    int i = _snp_name_map[snp];

    for (int p = 0; p < _fid.size(); p++)
    {
        int p_geno = 0;
        if (!_snp_1[i][p]) {
            p_geno |= (1 << 0);        
        } 
        if (!_snp_2[i][p]) {
            p_geno |= (1 << 1);
        }
        outfile << p_geno << endl;
    }
    outfile.close();
}

void gcta::dump_maternal_prob(string filename, string snp)
{
    ofstream outfile;
    outfile.open(filename);

    if (_snp_name_map.find(snp) == _snp_name_map.end()) {
        throw ("No such SNP found!");
    } 

    int col = _snp_name_map[snp];
    outfile << snp << endl;
    
    for (int row = 0; row < _keep.size(); row++)
    {
        outfile << _mat_geno(row, col) << endl; 
    }
    
    outfile.close();

}

#if 0 
float gcta::infer_maternal_genotype(string fid, string snp)
{
    if (_mu.empty())
    {
        calcu_mu();
        calcu_maf();
    }
    if (fid == "") 
        throw ("Need to supply a valid FID\n");

    //find the snp
    map<string, int>::iterator snp_iter = _snp_name_map.find(snp);
    if(snp_iter == _snp_name_map.end())
        throw ("Need to supply a valid SNP. SNP not found\n");

    parent_map_t parent_map; 

    // iterate over all individuals
    for (int i = 0; i < _fid.size(); i++)
    {
        string parent_id = _fa_id[i] + ":" + _mo_id[i];
        parent_map_t::iterator got = parent_map.find(parent_id);
        
        // match fid - same family 
        if (fid == _fid[i])
        {
            if (got == parent_map.end())
            {
                std::vector<int> index;
                index.push_back(i);
                parent_map[parent_id] = index;
            } else 
            {
                parent_map[parent_id].push_back(i);
            }
        }
    }

    int snp_index = snp_iter->second;
    float maf = _maf[snp_index];

    for (auto p = parent_map.begin(); p != parent_map.end(); p++)
    {
        std::vector<int> indices = p->second; 
        //means we have more than 2 individuals with the same parents
        if (indices.size() >= 2) 
        {
            //just take the first 2 individuals

            //person1 and person2 with same parent 
            int p1 = indices[0];
            int p2 = indices[1];

            int p1_geno = 0, p2_geno = 0;
            if (!_snp_1[snp_index][p1]) {
                p1_geno |= (1 << 0);        
            } 
            if (!_snp_2[snp_index][p1]) {
                p1_geno |= (1 << 1);
            }

            if (!_snp_1[snp_index][p2]) {
                p2_geno |= (1 << 0);        
            } 
            if (!_snp_2[snp_index][p2]) {
                p2_geno |= (1 << 1);
            }

            //cout << (enum SNPGenotype) p1_geno << " " << (enum SNPGenotype) p2_geno << endl;
            std::vector<float> result = calc_parental_prob((enum SNPGenotype) p1_geno, (enum SNPGenotype) p2_geno, 1- maf);
            //std::cout << result.size() << std::endl;
            //list_all_mat_probs(1-maf);
            float temp = result[1] + 2 * result[2];
            //printf("result: %f\n", temp);
            return result.size() ? result[1] + 2 * result[2] : 0.0f;
        }
    }

    //cout << "JUST A TEST" << endl;
    /*cout << "fid=" << fid << ", snp=" << snp << ", maf=" << maf << endl;  
    for (int i = 0; i < result.size(); i++) 
    {
        for (int j = 0; j < result[i].size(); j++)
        {
            cout << result[i][j] << ", ";
        }
        cout << endl;
    }*/
    return 0.0f;
}
#endif

/*eigenMatrix gcta::generate_maternal_probabilities()
{
    eigenMatrix result;
    result.resize(_keep.size(), _snp_num);
    int j;

    //todo: make this more efficient
    for (int i = 0; i < _keep.size(); i++)
    {
        #pragma omp parallel for private(i)
        for (j = 0; j < _snp_num; j++) 
        {
            string snp_name = _snp_name[j];
            string fid = _fid[i];
            result(i, j) = infer_maternal_genotype(fid, snp_name); 
            //cout << result(i, j) << ", ";
        }
        cout << i << "th completed" << endl;   
        //_snp_name_map.find(_snp_name[i])
    }
    return result;
}*/

void gcta::generate_maternal_probabilities()
{
    if (_mu.empty())
    {
        calcu_mu();
        calcu_maf();
    }
    parent_map_t parent_map; 

    // iterate over each individual
    for (int i = 0; i < _fid.size(); i++)
    {
        string parent_id = _fa_id[i] + ":" + _mo_id[i];
        parent_map_t::iterator got = parent_map.find(parent_id);
        
        if (got == parent_map.end())
        {
            // not found in parent_map, so we add vector  
            std::vector<int> index;
            index.push_back(i);
            parent_map[parent_id] = index;
        } else 
        {
            parent_map[parent_id].push_back(i);
        }
    }    

    _mat_geno.resize(_keep.size(), _snp_num);
    //int count1 = 0, count2 = 0;
    for (auto p = parent_map.begin(); p != parent_map.end(); p++)
    {
        std::vector<int> indices = p->second; 
        eigenMatrix probs;
        if (indices.size() >= 2) 
        {
            // we have more than 2 individuals with the same parents
            // take first two individuals
            int p1 = indices[0];
            int p2 = indices[1];
            calc_parental_prob_snps(p1, p2, probs);       
        } else {
            // person has no siblings   
            probs = eigenMatrix::Zero(1, _snp_num);
        }

        // update siblings with same probability
        for (int i = 0; i < indices.size(); i++)
        {
            _mat_geno.block(indices[i], 0, 1, _snp_num) = probs; 
        }
    }

    #if 0
    ofstream outfile;
    outfile.open("maternal_prob_dump.txt");
    //outfile << _mat_geno.block(0, 0, _fid.size(), 15);
    for (int i = 0; i < min(_snp_num, 30); i++) 
    {
        outfile << _snp_name[i];
        if (i != min(_snp_num, 30)-1) {
            outfile << ",";
        } 
    }
    outfile << endl;
    for (int row = 0; row < _keep.size(); row++)
    {
        for (int col = 0; col < min(_snp_num, 30); col++)
        {
            outfile << _mat_geno(row, col); 
            if (col != min(_snp_num, 30)-1) {
                outfile << ",";
            }
        }
        outfile << endl;
    }
    
    outfile.close();
    #endif
    //dump_genotypes("geno_dump.txt", "rs143500173");
    //dump_maternal_prob("maternal_prob_dump.txt",  "rs143500173");
}