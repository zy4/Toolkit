/*
 * Copyright 2018 Dept. of Protein Evolution, Max Planck Institute for Biology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package de.proteinevolution.tools.parameters

import de.proteinevolution.tel.TEL
import de.proteinevolution.tools.parameters.Parameter.{
  HHpredSelectsParameter,
  ModellerParameter,
  NumberParameter,
  SelectOption,
  SelectParameter,
  TextAreaParameter,
  TextInputParameter
}
import de.proteinevolution.tools.parameters.TextAreaInputType.TextAreaInputType
import javax.inject.{ Inject, Singleton }

/**
 * Provides the specification of the Parameters as they appear in the individual tools
 */
@Singleton
class ParamAccess @Inject() (tel: TEL) {

  def select(
      name: String,
      default: Option[String],
      maxSelectedOptions: Int = 1,
      onDetectedMSA: Option[String] = None
  ): SelectParameter =
    SelectParameter(
      name,
      default,
      tel.generateValues(name).toSeq.map(option => SelectOption(option._1, option._2)),
      maxSelectedOptions,
      onDetectedMSA
    )

  final val alignmentFormats = Seq(
    "fas" -> "fas",
    "a2m" -> "a2m",
    "a3m" -> "a3m",
    "sto" -> "sto",
    "psi" -> "psi",
    "clu" -> "clu"
  )

  final val samCCHelixRegex: Option[String] = Some("^[a-r];[a-zA-Z0-9];\\d+;\\d+$")

  def getParam(
      paramName: String,
      placeholderKey: String = "",
      sampleInputKey: String = "",
      alignmentInputType: TextAreaInputType = TextAreaInputType.SEQUENCE
  ): Parameter = paramName match {

    // Common
    case "ALIGNMENT"         => TextAreaParameter("alignment", alignmentInputType, placeholderKey, sampleInputKey)
    case "STANDARD_DB"       => select("standarddb", default = Some("alphafold_uniprot50"))
    case "PMIN"              => select("pmin", default = Some("20"))
    case "DESC"              => select("desc", default = Some("250"))
    case "MAXROUNDS"         => select("maxrounds", default = Some("1"))
    case "EVALUE"            => select("evalue", default = Some("1e-3"))
    case "OUTPUT_ORDER"      => select("output_order", default = Some("input"))
    case "PCOILS_INPUT_MODE" => select("pcoils_input_mode", default = Some("0"), onDetectedMSA = Some("1"))
    case "MATRIX"            => select("matrix", default = Some("BLOSUM62"))
    case "MIN_SEQID_QUERY"   => select("min_seqid_query", default = Some("0"))

    // DIAMOND-DeepClust
    case "DIAMOND_MIN_SEQID"   => select("diamond_min_seqid", default = Some("80"))
    case "DIAMOND_MIN_ALN_COV" => select("diamond_min_aln_cov", default = Some("80"))

    // HHblits
    case "HHBLITSDB" => select("hhblitsdb", default = Some("UniRef30"))
    case "HHBLITS_INCL_EVAL" =>
      select("hhblits_incl_eval", default = Some("1e-3"))

    // pLM-BLAST
    case "PLMBLASTDB"               => select("plmblastdb", default = Some("ECOD30"))
    case "COSINE_PERCENTILE_CUTOFF" => select("cosine_percentile_cutoff", default = Some("70"))
    case "ALIGNMENT_CUTOFF"         => select("alignment_cutoff", default = Some("0.30"))
    case "WIN_LEN"                  => select("win_len", default = Some("15"))
    case "MERGE_HITS" => select("merge_hits", default = Some("1"))
    case "PLM_ALN_MODE" => select("plm_aln_mode", default = Some("loc"))
    case "SPAN" => select("span", default = Some("25"))
    case "SIGMA_FACTOR" => select("sigma_factor", default = Some("2"))
    case "BFACTOR" => select("bfactor", default = Some("2"))
    case "PLM_GAP_EXT" => select("plm_gap_ext", default = Some("0.5"))



    // HHpred
    case "TWOTEXTALIGNMENT" =>
      TextAreaParameter("alignment", alignmentInputType, placeholderKey, sampleInputKey, allowsTwoTextAreas = true)
    case "HHPRED_DB_PROTEOMES" =>
      HHpredSelectsParameter(
        "hhsuitedb",
        tel.generateValues("hhsuitedb").toSeq.map(option => SelectOption(option._1, option._2)),
        "proteomes",
        tel.generateValues("proteomes").toSeq.map(option => SelectOption(option._1, option._2)),
        maxSelectedOptions = 4,
        default = Some("mmcif70/pdb70"),
        defaultProteomes = None
      )
    case "MSA_GEN_METHOD"   => select("msa_gen_method", default = Some("UniRef30"))
    case "MSA_GEN_MAX_ITER" => select("msa_gen_max_iter", default = Some("3"), onDetectedMSA = Some("0"))
    case "HHPRED_INCL_EVAL" =>
      select("hhpred_incl_eval", default = Some("1e-3"))
    case "MIN_COV"      => select("min_cov", default = Some("20"))
    case "SS_SCORING"   => select("ss_scoring", default = Some("2"))
    case "ALIGNMACMODE" => select("alignmacmode", default = Some("loc"))
    case "MACTHRESHOLD" => select("macthreshold", default = Some("0.3"))

    // HMMER
    case "HMMER_DB" => select("hmmerdb", default = Some("alphafold_uniprot50"))
    case "MAX_HHBLITS_ITER" =>
      select("max_hhblits_iter", default = Some("1"), onDetectedMSA = Some("0"))

    // PatternSearch
    case "PATSEARCH_DB" => select("patsearchdb", default = Some("alphafold_uniprot50"))
    case "GRAMMAR"      => select("grammar", default = Some("pro"))
    case "SEQCOUNT"     => select("seqcount", default = Some("500"))

    // PSI-BLAST
    case "BLAST_INCL_EVAL" =>
      select("blast_incl_eval", default = Some("1e-3"))

    // Kalign
    case "GAP_OPEN" => NumberParameter("gap_open", default = Some(11))
    case "GAP_EXT_KALN" =>
      NumberParameter(
        "gap_ext_kaln",
        min = Some(0),
        max = Some(10),
        step = Some(0.01),
        default = Some(0.85)
      )
    case "GAP_TERM" =>
      NumberParameter(
        "gap_term",
        min = Some(0),
        max = Some(10),
        step = Some(0.01),
        default = Some(0.45)
      )
    case "BONUSSCORE" =>
      NumberParameter("bonusscore", min = Some(0), max = Some(10), step = Some(0.01), default = Some(0))

    // MAFFT
    case "MAFFT_GAP_OPEN" =>
      NumberParameter(
        "mafft_gap_open",
        min = Some(0),
        max = Some(10),
        step = Some(0.01),
        default = Some(1.53)
      )
    case "OFFSET" =>
      NumberParameter("offset", min = Some(0), max = Some(10), step = Some(0.01), default = Some(0.0))

    // HHrepID
    case "MSA_GEN_MAX_ITER_HHREPID" =>
      select("msa_gen_max_iter_hhrepid", default = Some("3"), onDetectedMSA = Some("0"))
    case "SCORE_SS"           => select("score_ss", default = Some("2"))
    case "REP_PVAL_THRESHOLD" => select("rep_pval_threshold", default = Some("1e-2"))
    case "SELF_ALN_PVAL_THRESHOLD" =>
      select("self_aln_pval_threshold", default = Some("1e-1"))
    case "MERGE_ITERS"            => select("merge_iters", default = Some("3"))
    case "DOMAIN_BOUND_DETECTION" => select("domain_bound_detection", default = Some("1"))

    // MARCOIL
    case "MATRIX_MARCOIL"         => select("matrix_marcoil", default = Some("mtk"))
    case "TRANSITION_PROBABILITY" => select("transition_probability", default = Some("1"))

    // PCOILS
    case "PCOILS_WEIGHTING" => select("pcoils_weighting", default = Some("1"))
    case "PCOILS_MATRIX"    => select("pcoils_matrix", default = Some("2"))

    // REPPER
    case "REPPER_INPUT_MODE" => select("repper_input_mode", default = Some("0"))
    case "WINDOW_SIZE"       => NumberParameter("window_size", default = Some(100))
    case "PERIODICITY_MIN"   => NumberParameter("periodicity_min", default = Some(2))
    case "PERIODICITY_MAX"   => NumberParameter("periodicity_max", default = Some(100))
    case "FTWIN_THRESHOLD"   => NumberParameter("ftwin_threshold", default = Some(6))
    case "REPWIN_THRESHOLD"  => NumberParameter("repwin_threshold", default = Some(2))

    // TPRpred
    case "EVAL_TPR" => select("eval_tpr", default = Some("1e-2"))

    // Ali2D
    case "INVOKE_PSIPRED" =>
      NumberParameter(
        "invoke_psipred",
        min = Some(0),
        max = Some(100),
        default = Some(30)
      )

    // HHomp
    case "HHOMPDB"   => select("hhompdb", default = Some("HHompDB_1.0.hhm"))
    case "ALIGNMODE" => select("alignmode", default = Some("loc"))

    // Quick2D
    case "TARGET_PSI_DB" => select("target_psi_db", default = Some("nr90"))
    case "QUICK_ITERS"   => select("quick_iters", default = Some("3"))

    // MODELLER
    case "REGKEY" => ModellerParameter("regkey", "Enter MODELLER-key (see help pages for details)")

    // SamCC
    case "SAMCC_HELIXONE" =>
      TextInputParameter(
        "samcc_helixone",
        "CC_first_position;chain;start_pos;end_pos",
        samCCHelixRegex,
        Some("a;A;2;30"),
        disableRemember = Some(true)
      )
    case "SAMCC_HELIXTWO" =>
      TextInputParameter(
        "samcc_helixtwo",
        "CC_first_position;chain;start_pos;end_pos",
        samCCHelixRegex,
        Some("a;B;2;30"),
        disableRemember = Some(true)
      )
    case "SAMCC_HELIXTHREE" =>
      TextInputParameter(
        "samcc_helixthree",
        "CC_first_position;chain;start_pos;end_pos",
        samCCHelixRegex,
        Some("a;C;2;30"),
        disableRemember = Some(true)
      )
    case "SAMCC_HELIXFOUR" =>
      TextInputParameter(
        "samcc_helixfour",
        "CC_first_position;chain;start_pos;end_pos",
        samCCHelixRegex,
        Some("a;D;2;30"),
        disableRemember = Some(true)
      )
    case "SAMCC_PERIODICITY" => select("samcc_periodicity", default = Some("7"))
    case "EFF_CRICK_ANGLE"   => select("eff_crick_angle", default = Some("1"))

    // CLANS
    case "CLANS_EVAL" => select("clans_eval", default = Some("1e-4"))

    // MMseqs2
    case "MIN_SEQID"       => select("min_seqid", default = Some("0.8"))
    case "MIN_ALN_COV"     => select("min_aln_cov", default = Some("0.8"))
    case "CLUSTERING_MODE" => select("clustering_mode", default = Some("cluster"))

    // PhyML
    case "MATRIX_PHYML"  => select("matrix_phyml", default = Some("LG"))
    case "NO_REPLICATES" => select("no_replicates", default = Some("0"))

    // 6FrameTranslation
    case "INC_NUCL"       => select("inc_nucl", default = Some("t"))
    case "AMINO_NUCL_REL" => select("amino_nucl_rel", default = Some("t"))
    case "CODON_TABLE"    => select("codon_table", default = Some("1"))

    // BackTranslator
    case "INC_AMINO"            => select("inc_amino", default = Some("fas"))
    case "GENETIC_CODE"         => select("genetic_code", default = Some("1"))
    case "CODON_TABLE_ORGANISM" => TextInputParameter("codon_table_organism", inputPlaceholder = "")

    // FormatSeq
    case "IN_FORMAT"  => select("in_format", default = Some("a3m"))
    case "OUT_FORMAT" => select("out_format", default = Some("fas"))

    // HHfilter
    case "MAX_SEQID" => NumberParameter("max_seqid", default = Some(90))
    case "MIN_QUERY_COV" =>
      NumberParameter(
        "min_query_cov",
        min = Some(0),
        max = Some(100),
        default = Some(0)
      )
    case "NUM_SEQS_EXTRACT" =>
      NumberParameter("num_seqs_extract", default = Some(0))
  }

  val paramGroups: Map[String, Seq[String]] = Map(
    "Input" -> Seq(
      getParam("ALIGNMENT").name,
      getParam("STANDARD_DB").name,
      getParam("HHPRED_DB_PROTEOMES").name,
      getParam("HHBLITSDB").name,
      getParam("PLMBLASTDB").name,
      getParam("HHOMPDB").name,
      getParam("HMMER_DB").name,
      getParam("PATSEARCH_DB").name,
      getParam("REGKEY").name,
      getParam("GRAMMAR").name,
      getParam("SAMCC_HELIXONE").name,
      getParam("SAMCC_HELIXTWO").name,
      getParam("SAMCC_HELIXTHREE").name,
      getParam("SAMCC_HELIXFOUR").name,
      getParam("TARGET_PSI_DB").name,
      getParam("QUICK_ITERS").name,
      getParam("PCOILS_INPUT_MODE").name,
      getParam("REPPER_INPUT_MODE").name,
      getParam("IN_FORMAT").name,
      getParam("OUT_FORMAT").name
    )
  )
}
