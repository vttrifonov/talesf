from libc.stdlib cimport malloc, free, calloc
from libc.stdio cimport FILE, stdout, fopen, fclose
import re

ctypedef void (*valuefreefunc)(void *)

cdef extern from "Hashmap.h":

    ctypedef struct Hashmap:
        pass

    Hashmap *hashmap_new(int size)
    int hashmap_add(Hashmap *hash, char *key, void *value)
    void hashmap_delete(Hashmap *hash, valuefreefunc)
    int hashmap_size(Hashmap *hash)
    void *hashmap_get(Hashmap *hash, char *key)
    char **hashmap_keys(Hashmap *hash)

cdef extern from "Array.h":
        
    ctypedef struct Array:
        pass
    
    Array *array_new()
    void array_add(Array *r, void *e)
    void *array_get(Array *r, int index)
    void array_delete(Array *r, valuefreefunc)

cdef extern from "bcutils.h":
    double *double_array(double a, double c, double g, double t, double dummy)
    Hashmap *get_diresidue_probabilities(Array *rvdseq, double w)
    Hashmap *convert_probabilities_to_scores(Hashmap *diresidue_probabilities)

cdef get_best_score(rvd_seq, Hashmap *rvdscores):
    
    cdef:
        int i,j
        double best_score = 0.0
        double min_score = -1.0
        double *scores
    
    for i in range(len(rvd_seq)):
        scores = <double*> hashmap_get(rvdscores, rvd_seq[i])
        if scores == NULL:
            return -1.0
        for j in range(4):
            if j == 0 or scores[j] < min_score:
                min_score = scores[j]
        best_score += min_score
    
    return best_score

cdef extern from "talesf.h":
    int run_talesf_task(Hashmap *kwargs)

def ScoreTalesfTask(char *seqfilename, rvd_string, char *output_filepath, char *log_filepath, bint forwardonly, int c_upstream, double cutoff, int numprocs, char *organism_name):
    
    cdef:
        int i, j
        double weight = 0.9
        
        Hashmap *talesf_kwargs = hashmap_new(32)
        
        int BIGGEST_RVD_SCORE_EVER = 100
        
        Array *rvd_array = array_new()
        
        char *rvd_str
        
        char *rvd_string_ptr = rvd_string
    
    split_rvd_string = re.split(' |_', rvd_string)
    
    for rvd in split_rvd_string:
        rvd_str = rvd
        array_add(rvd_array, rvd_str)
    
    cdef:
        
        Hashmap *diresidue_probabilities = get_diresidue_probabilities(rvd_array, weight)
        Hashmap *diresidue_scores = convert_probabilities_to_scores(diresidue_probabilities)
    
    hashmap_delete(diresidue_probabilities, NULL)
    hashmap_add(diresidue_scores, "XX", double_array(0, 0, 0, 0, BIGGEST_RVD_SCORE_EVER))
    
    cdef:
        double **scoring_matrix = <double**> calloc(hashmap_size(diresidue_scores), sizeof(double*))
        
        unsigned int *rvd_seq = <unsigned int*> calloc(len(split_rvd_string), sizeof(unsigned int))
        
        unsigned int rvd_seq_len = len(split_rvd_string)
        
        double best_score = get_best_score(split_rvd_string, diresidue_scores)
    
    cdef char **diresidues = hashmap_keys(diresidue_scores)
    
    rvd_to_int = {}
    
    for i in range(hashmap_size(diresidue_scores)):
        rvd_to_int[diresidues[i]] = i
        scoring_matrix[i] = <double*> hashmap_get(diresidue_scores, diresidues[i])
        scoring_matrix[i][4] = BIGGEST_RVD_SCORE_EVER
    
    for i in range(rvd_seq_len):
        rvd_seq[i] = rvd_to_int[split_rvd_string[i]]
    
    hashmap_add(talesf_kwargs, "seq_filename", seqfilename)
    hashmap_add(talesf_kwargs, "rvd_seq", rvd_seq)
    hashmap_add(talesf_kwargs, "rvd_seq_len", &rvd_seq_len)
    hashmap_add(talesf_kwargs, "rvd_string", rvd_string_ptr)
    hashmap_add(talesf_kwargs, "best_score", &best_score)
    hashmap_add(talesf_kwargs, "scoring_matrix", scoring_matrix)
    hashmap_add(talesf_kwargs, "output_filepath", output_filepath)
    hashmap_add(talesf_kwargs, "log_filepath", log_filepath)
    hashmap_add(talesf_kwargs, "weight", &weight)
    hashmap_add(talesf_kwargs, "cutoff", &cutoff)
    hashmap_add(talesf_kwargs, "c_upstream", &c_upstream)
    hashmap_add(talesf_kwargs, "num_procs", &numprocs)
    hashmap_add(talesf_kwargs, "organism_name", organism_name)
    
    hashmap_add(talesf_kwargs, "forward_only", &forwardonly)
    
    cdef int task_result = run_talesf_task(talesf_kwargs)
    
    free(scoring_matrix)
    free(diresidues)
    array_delete(rvd_array, NULL)
    hashmap_delete(diresidue_scores, free)
    free(rvd_seq)
    
    hashmap_delete(talesf_kwargs, NULL)
    
    return task_result
