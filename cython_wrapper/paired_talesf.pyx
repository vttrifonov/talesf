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

cdef extern from "pairedtalesf.h":
    int run_paired_talesf_task(Hashmap *kwargs)

def ScorePairedTalesfTask(char *seqfilename, rvd_string, rvd_string2, char *output_filepath, char *log_filepath, int c_upstream, double cutoff, int spacer_min, int spacer_max, int numprocs, char *organism_name):
    
    cdef:
        int i, j
        double weight = 0.9
        int count_only = 0
        
        Hashmap *paired_talesf_kwargs = hashmap_new(32)
        
        int BIGGEST_RVD_SCORE_EVER = 100
        
        Array *rvd_array = array_new()
        
        char *rvd_str
        
        char *rvd_string_ptr = rvd_string
        char *rvd_string2_ptr = rvd_string2
    
    split_rvd_string = re.split(' |_', rvd_string)
    split_rvd_string2 = re.split(' |_', rvd_string2)
    
    joint_split_rvd_strings = split_rvd_string + split_rvd_string2
    
    for rvd in joint_split_rvd_strings:
        rvd_str = rvd
        array_add(rvd_array, rvd_str)
    
    cdef:
        
        Hashmap *diresidue_probabilities = get_diresidue_probabilities(rvd_array, weight)
        Hashmap *diresidue_scores = convert_probabilities_to_scores(diresidue_probabilities)
        double **scoring_matrix = <double**> calloc(hashmap_size(diresidue_scores), sizeof(double*))
        
        unsigned int *rvd_seq = <unsigned int*> calloc(len(split_rvd_string), sizeof(unsigned int))
        unsigned int *rvd_seq2 = <unsigned int*> calloc(len(split_rvd_string2), sizeof(unsigned int))
        unsigned int **rvd_seqs = [rvd_seq, rvd_seq2]
        
        unsigned int rvd_seq_len = len(split_rvd_string)
        unsigned int rvd_seq2_len = len(split_rvd_string2)
        unsigned int *rvd_seqs_lens = [rvd_seq_len, rvd_seq2_len]
        
        double best_score = get_best_score(split_rvd_string, diresidue_scores)
        double best_score2 = get_best_score(split_rvd_string2, diresidue_scores)
        double *best_scores = [best_score, best_score2]
    
    hashmap_delete(diresidue_probabilities, NULL)
    hashmap_add(diresidue_scores, "XX", double_array(0, 0, 0, 0, BIGGEST_RVD_SCORE_EVER))
    cdef char **diresidues = hashmap_keys(diresidue_scores)
    
    rvd_to_int = {}
    
    for i in range(hashmap_size(diresidue_scores)):
        rvd_to_int[diresidues[i]] = i
        scoring_matrix[i] = <double*> hashmap_get(diresidue_scores, diresidues[i])
        scoring_matrix[i][4] = BIGGEST_RVD_SCORE_EVER
    
    for i in range(rvd_seq_len):
        rvd_seq[i] = rvd_to_int[split_rvd_string[i]]
    
    for i in range(rvd_seq2_len):
        rvd_seq2[i] = rvd_to_int[split_rvd_string2[i]]
    
    hashmap_add(paired_talesf_kwargs, "seq_filename", seqfilename)
    hashmap_add(paired_talesf_kwargs, "rvd_seqs", rvd_seqs)
    hashmap_add(paired_talesf_kwargs, "rvd_seqs_lens", rvd_seqs_lens)
    hashmap_add(paired_talesf_kwargs, "rvd_string", rvd_string_ptr)
    hashmap_add(paired_talesf_kwargs, "rvd_string2", rvd_string2_ptr)
    hashmap_add(paired_talesf_kwargs, "best_scores", best_scores)
    hashmap_add(paired_talesf_kwargs, "scoring_matrix", scoring_matrix)
    hashmap_add(paired_talesf_kwargs, "output_filepath", output_filepath)
    hashmap_add(paired_talesf_kwargs, "log_filepath", log_filepath)
    hashmap_add(paired_talesf_kwargs, "weight", &weight)
    hashmap_add(paired_talesf_kwargs, "cutoff", &cutoff)
    hashmap_add(paired_talesf_kwargs, "c_upstream", &c_upstream)
    hashmap_add(paired_talesf_kwargs, "spacer_min", &spacer_min)
    hashmap_add(paired_talesf_kwargs, "spacer_max", &spacer_max)
    hashmap_add(paired_talesf_kwargs, "num_procs", &numprocs)
    hashmap_add(paired_talesf_kwargs, "organism_name", organism_name)
    
    hashmap_add(paired_talesf_kwargs, "count_only", &count_only)
    
    cdef int task_result = run_paired_talesf_task(paired_talesf_kwargs)
    
    free(scoring_matrix)
    free(diresidues)
    array_delete(rvd_array, NULL)
    hashmap_delete(diresidue_scores, free)
    free(rvd_seq)
    free(rvd_seq2)
    
    
    return task_result
