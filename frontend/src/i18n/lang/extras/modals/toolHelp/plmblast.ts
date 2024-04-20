/* tslint:disable:max-line-length */

export default {
    en: {
        toolHelpModals: {
            plmblast: {
                overview: `
                <p>pLM-BLAST is a sensitive remote homology detection tool that is based on the comparison of sequence 
                embeddings obtained from the protein language model ProtT5. It is available as an easy-to-use web server 
                within the MPI Bioinformatics Toolkit but can be also used as a standalone package, which is available 
                at <a href = https://github.com/labstructbioinf/pLM-BLAST target="_blank" rel="noopener">
                https://github.com/labstructbioinf/pLM-BLAST. </a> Pre-computer databases can be downloaded from 
                <a href = http://ftp.tuebingen.mpg.de/pub/protevo/toolkit/databases/plmblast_dbs target="_blank" rel="noopener">
                 http://ftp.tuebingen.mpg.de/pub/protevo/toolkit/databases/plmblast_dbs</a>.</p>`,
                parameters: [
                    {
                        title: 'Input',
                        content: `You have the option to enter a single protein sequence FASTA.<br/><br/>
                @:toolHelpModals.common.singleseq`,
                    },
                    {
                        title: 'Select target database',
                        content: `<p>Select domain database(s) of template embeddings against which you want to compare the query.</p>
                    <em>ECOD30, ECOD70</em>
                    <p> These databases are versions of the <a href = http://prodata.swmed.edu/ecod/ target="_blank" rel="noopener">
                        Evolutionary Classification of Protein Domains (ECOD) database</a> filtered for a maximum of 30%,
                        50%, and 70% sequence identity, respectively.</p>`,
                    },
                    {
                        title: 'Cosine similarity percentile cut-off',
                        content: `A pre-screening procedure can be used to improve the performance of the database search.
                        First, the database of flattened embeddings is searched using the flattened query embedding and the
                        cosine similarity metric, and then the actual pLM-BLAST comparisons are performed only for matches
                        above the user-specified cut-off. The cut-off is expressed as the n-th percentile of all cosine
                        similarity scores. The higher the pre-screening cut-off, the faster and less sensitive the search
                        will be, and vice versa.`,
                    },
                    {
                        title: 'Alignment score cut-off',
                        content: `Each query-target alignment is assigned a score from 0 to 1. The alignment cut-off defines
                        the minimum score for reporting a match. The higher the cut-off, the more stringent the search. Also
                        note that only matches that have passed the pre-filtering step (see "Cosine Similarity Percentile Cutoff")
                        are considered.`,
                    },
                    {
                        title: 'Window length',
                        content: `A moving average is used to detect local alignments within the paths determined by the traceback
                        procedure. Increasing the window size results in longer local alignments, but may decrease sensitivity.
                        This parameter is not used in Global Alignment mode.`,
                    },
                    {
                        title: 'Merge hits',
                        content: `Since pLM-BLAST tends to return rather short alignments, an optional procedure can be 
                        used to merge matches to a single database entry. Such a merged match can consist of two or more 
                        unmerged hits, and its score, similarity, and identity are defined as the average of the values 
                        from the individual sub-hits.`,
                    },
                    {
                        title: 'Max target hits',
                        content: `This parameter controls how many matches are displayed in the results. This filter is
                        applied after the alignment score cutoff filter.`,
                    },
                    {
                        title: 'Alignment mode',
                        content: `Specifies whether to return local or global alignments.`,
                    },
                    {
                        title: 'Minimal hit span',
                        content: `Specifies the minimum length of matches returned. The minimum match length cannot be
                        less than the window length.`,
                    },
                    {
                        title: 'Sigma factor',
                        content: `The sigma factor defines the cutoff at which the background signal is discarded when 
                        searching for significant local alignments. Increasing (>2) or decreasing (<2) the cutoff makes 
                        the algorithm stricter or more permissive, respectively. This parameter is not used in Global 
                        Alignment mode.`,
                    },
                    {
                        title: 'Bfactor',
                        content: `Using a Bfactor value greater than 1 will reduce the initial number of alignments by
                        ignoring alignments that are very close to each other, thus increasing the search speed. This
                        parameter is not used in Global Alignment mode.`,
		    },
                    {
                        title: 'Gap extension penalty',
                        content: `Setting this parameter to a value greater than 0 penalizes the introduction of gaps
			in the generated alignments.`,
	            },

                ],
                references: `<p>pLM-BLAST – distant homology detection based on direct comparison of sequence 
                                representations from protein language models.<br>
                                Kaminski K, Ludwiczak K, Alva V, Dunin-Horkawicz S. bioRxiv. 2022.</p>`,
            },
        },
    },
};
