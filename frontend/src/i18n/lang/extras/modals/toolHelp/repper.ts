/* tslint:disable:max-line-length */

export default {
    en: {
        toolHelpModals: {
            repper: {
                overview: `<p>REPPER is a server that detects and analyses regions with short gapless REPeats in protein sequences.
            It finds PERiodicities by Fourier Transform (FTwin) and internal homology analysis (REPwin). The exact
            locations of the periodic patterns are traced by a sliding window. The output is complemented by coiled coil
            prediction (COILS) and optionally by secondary structure prediction (PSIPRED).</p>
        <h5><b>FTwin</b></h5>
        <p>Sequence analysis by a sliding window Fourier transform (FT) method calculating the periodicity spectra for
            all possible sequence windows. The overview graph shows all significant periodicities and where in the
            sequence they are located. Periodicities are considered to be significant if their FT intensity is above
            the given threshold.</p>
        <p>Clicking on a colored line on the output graph provides the complete FT spectrum for the sequence fragment
            defined by the line (the fragment may include closely neighbored lines as well). Alternatively, an arbitrary
            sequence range may be given (in the boxes under the output graph) in order to get a FT spectrum of this
            particular sequence fragment. <a href = http://ftp.tuebingen.mpg.de/ebio/protevo/FTwin target="_blank" rel="noopener">Download</a></p>
        <h5><b>REPwin</b></h5>
        <p>REPWin compares a protein sequence with itself within a sliding window of choosable length W. For each
            window starting position i and periodicity p it calculates</p>
        <div align = center><em>Score(i,p) = Sum_k Sum_j S(x[j],x[j+kp])</em></div><br/>
        <p>(x[j],x[j+kp]) is the Gonnet substitution matrix element for residues x[j] and x[j+kp]. The sum runs over
            all k and j such that j and j+kp are inside the window [i,..,i+W-1]. The score S(i,p) is normalized by
            dividing through the standard deviation of S(i,p) for nonperiodic sequences. The final score value for each
            residue number i and periodicity p is the maximum over all windows containing residue i.</p>
        <p>Clicking on a colored line on the output graph displays the corresponding alignment that shows the internal
            homology. <a href = "http://ftp.tuebingen.mpg.de/ebio/protevo/REPwin/" target="_blank" rel="noopener">Download</a></p>
        <h5><b>COILS</b></h5>
        <p>COILS is a program that compares a sequence to a database of known parallel two-stranded coiled-coils and
            derives a similarity score. By comparing this score to the distribution of scores in globular and
            coiled-coil proteins, the program then calculates the probability that the sequence will adopt a
            coiled-coil conformation. COILS is described in <a href = https://www.ncbi.nlm.nih.gov/pubmed/2031185?dopt=Abstract target="_blank" rel="noopener">Lupas, Science 252: 1162-1164 (1991)</a> and Lupas, Meth. Enzymology 266: 513-525 (1996). It is based on a prediction protocol proposed by <a href = https://www.ncbi.nlm.nih.gov/pubmed/7165792?dopt=Abstract>Parry, Biosci. Rep. 2: 1017-1024 (1982)</a>. You can also <a href = http://ftp.tuebingen.mpg.de/ebio/protevo/COILS target="_blank" rel="noopener">download</a> the software.</p>
        <h5><b>PSIPRED</b></h5>
        <p>PSIPRED predicts the secondary structure as described in <a href = https://www.ncbi.nlm.nih.gov/pubmed/10493868?dopt=Abstract target="_blank" rel="noopener">Jones, J.Mol.Biol. 292:195-202 (1999)</a>.</p>`,
                parameters: [
                    {
                        title: 'Input',
                        content: `You have the option to enter or upload either a single protein sequence or a multiple sequence
                    alignment (MSA) in FASTA or CLUSTAL format.<br/><br/>
                @:toolHelpModals.common.singleseq
                @:toolHelpModals.common.msa`,
                    },
                    {
                        title: 'Input mode',
                        content: `<p>Use input as given: Performs all calculations on the basis of a given sequence or a given
                    sequence alignment in FASTA format. Either put in a sequence or sequence alignment of your choice.</p>
                    <p>Run PSI-BLAST: Performs a PSI-BLAST run on the given sequence (two iterations, E-value for
                    inclusion 0.0001). Hits are aligned and filtered to a minimum sequence identity of 40% with
                    the query sequence and a minimum coverage of 20%. This alignment is then the basis for the
                    input profile. In case you are not happy with the used PSI-BLAST parameters, run PSI-BLAST
                    separately with the parameters of your choice, create a multiple alignment in FASTA format,
                    and use this as REPPER input.</p>`,
                    },
                    {
                        title: 'Window size (FTwin and REPwin)',
                        content: `Size of the sliding window. Must be smaller than sequence length.`,
                    },
                    {
                        title: 'Periodicity range (FTwin)',
                        content: `Periodicity range for Fourier Transform. Minimal periodicity must be at least 2,
            maximal periodicity should not exceed the window size (otherwise it is set to window size automatically).`,
                    },
                    {
                        title: 'Scale (FTwin)',
                        content: `A Fourier Transform cannot deal with amino acids. Amino acids are therefore represented
            by numbers. Depending on what one is looking for (e.g. repeat patterns of hydrophobic residues),
            those numbers represent a certain property (in this case hydrophobicity). The default values are
            set to the Kyte-Doolittle hydrophobicity scale, but you can also switch to an empirical set of
            binary hydrophobic values (standard binary) by clicking the button or you can manually put in values
            for your own purposes , either by changing each value manually or by pasting all amino acids into
            the prepared form (own scale).`,
                    },
                    {
                        title: 'FTwin threshold',
                        content: `<p>Refinement option. For each window, the FTwin output graph shows all periodicities above</p>
                    <p>average+threshold*standard_deviation,</p>
                    <p>while "average" is the average intensity of the W Fourier transform periodicities of the given
                    window (W = window size), "threshold" is the parameter given by the user, and "standard_deviation"
                    is the standard deviation of the intensities of the W Fourier transform periodicities of the given window.</p>
                    <p>Whether the threshold parameter is reasonable depends on the values that are attributed
                    to the amino acids, but also on the window size. Decrease this value if you see too much in the
                    FTwin graph, increase it if you see too little or nothing.</p>`,
                    },
                    {
                        title: 'REPwin threshold',
                        content: `Refinement option. The default value of 2 is reasonable for most applications.
             Decrease this value if you see too much in the REPwin graph, increase it if you see too little or nothing.`,
                    },
                ],
                references: `<p>Gruber M., Söding J., Lupas AN. (2005)
            <b>REPPER -- Repeats and their periodicities in fibrous proteins. </b>
            Nucl Acids Res 33(Web Server issue): W239-W243.
            <a href = https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1160166 target="_blank" rel="noopener">PMCID: PMC1160166</a></p>`,
            },
        },
    },
};
