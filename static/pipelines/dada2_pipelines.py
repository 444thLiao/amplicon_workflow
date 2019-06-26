from qiime2.plugins.dada2.methods import denoise_paired


def dada2_pipelines(dada2_input=None,
                    trunc_len_f=140,
                    trunc_len_r=150,
                    n_threads=0,  # all threads
                    trunc_q=2,  # default
                    n_reads_learn=1000000,  # default
                    max_ee=2.0,  # default
                    **kwargs):
    tab, rep, stats = denoise_paired(dada2_input,
                                     trunc_len_f=trunc_len_f,
                                     trunc_len_r=trunc_len_r,
                                     n_threads=n_threads,
                                     trunc_q=trunc_q,
                                     n_reads_learn=n_reads_learn,
                                     max_ee=max_ee
                                     )

    return (tab,
            rep,
            stats)
