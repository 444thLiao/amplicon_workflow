from qiime2.plugins.deblur.methods import denoise_16S


def deblur_pipelines(deblur_input=None,
                     trim_length=250,
                     sample_stats=True,
                     mean_error=0.005,  # default
                     indel_prob=0.01,  # default
                     indel_max=3,  # default
                     min_reads=10,  # default
                     min_size=2,  # default
                     jobs_to_start=7,
                     hashed_feature_ids=True,  # default
                     **kwargs):
    tab, rep, stats = denoise_16S(demultiplexed_seqs=deblur_input,
                                  trim_length=trim_length,
                                  sample_stats=sample_stats,
                                  mean_error=mean_error,  # default
                                  indel_prob=indel_prob,  # default
                                  indel_max=indel_max,  # default
                                  min_reads=min_reads,  # default
                                  min_size=min_size,  # default
                                  jobs_to_start=jobs_to_start,
                                  hashed_feature_ids=hashed_feature_ids,  # default
                                  )

    return (tab,
            rep,
            stats)
