export interface SearchAlignmentItemRender {
    query: any;
    agree: string;
    template: any;
}

export interface SearchAlignmentItem extends SearchAlignmentItemRender {
    num: number;
    acc: string;
    name: string;
}

export interface SearchHitItem {
    num: number;
}

export interface SearchHitsResponse {
    hits: SearchHitItem[];
    total: number;
    totalNoFilter: number;
    start: number;
    end: number;
}

export interface SearchAlignmentsResponse<T extends SearchAlignmentItem, S extends HHInfoResult> {
    alignments: T[];
    info: S;
    total: number;
    start: number;
    end: number;
}

export interface HHompAlignmentItem extends SearchAlignmentItem {
    alignedCols: number;
    probabHit: number;
    probabOMP: number;
    eval: number;
    score: number;
    identities: number;
}

export interface HHblitsAlignmentItem extends SearchAlignmentItem {
    alignedCols: number;
    probab: number;
    eval: number;
    score: number;
    identities: number;
}

export interface HHpredAlignmentItem extends SearchAlignmentItem {
    structLink: boolean;
    dbLink: string;
    probab: number;
    eval: number;
    score: number;
    ident: number;
    similarity: number;
    ssScore: number;
    alignedCols: number;
}

export interface HMMERAlignmentItem extends SearchAlignmentItem {
    fastaLink: string;
    fullEval: number;
    eval: number;
    bitScore: number;
    hitLen: number;
    observedDomains: number;
}

export interface PLMBLASTAlignmentItem extends SearchAlignmentItem {
    fastaLink: string;
    fullEval: number;
    eval: number;
    bitScore: number;
    hitLen: number;
    matchLen: number;
}

export interface PSIBLASTAlignmentItem extends SearchAlignmentItem {
    fastaLink: string;
    eval: number;
    bitScore: number;
    score: number;
    ident: number;
    perIdent: number;
    pos: number;
    perPos: number;
    gap: number;
    perGap: number;
    refLen: number;
    hitLen: number;
    fullEval: number;
}

export interface AlignmentResultResponse {
    alignments: AlignmentItem[];
    total: number;
    start: number;
    end: number;
}

export interface AlignmentItem {
    num: number;
    accession: string;
    seq: string;
}

export interface QueryItem {
    header: string;
    sequence: string;
}

export interface Quick2dResults {
    jobID: string;
    query: QueryItem;
    results: { [key: string]: string };
}

export interface PatsearchResults {
    jobID: string;
    results: {
        hits: PatsearchHit[];
        regex: string;
    };
}

export interface PatsearchHit {
    name: string;
    seq: string;
    matches: PatsearchMatch[];
}

export interface PatsearchMatch {
    i: number; // i: start index
    n: number; // n: length of match
}

export interface TprpredResults {
    desc: string[];
    hits: string[];
}

export interface HhrepidResults {
    jobID: string;
    results: {
        reptypes: HhrepidReptypes[];
    };
}

export interface HhrepidReptypes {
    pval: string;
    reps: HhrepidReptype[];
    len: number;
    typ: string;
    num: number;
}

export interface HhrepidReptype {
    prob: string;
    pval: string;
    loc: string;
    seq: string;
}

export interface ProbEvalList {
    type: string;
    vals: number[];
}

export interface HitMapResponse {
    hitAreas: HitMapItem[];
    queryLength: number;
    resubmitStart: number;
    resubmitEnd: number;
}

export interface HitMapItem {
    num: number;
    title: string;
    b: number;
    t: number;
    l: number;
    r: number;
}

export interface HHInfoResult {
    num_hits: number;
    query: {
        accession: string;
        seq: string;
    };
}

export interface PsiblastHHInfoResult extends HHInfoResult {
    tm: number;
    coil: number;
    signal: string;
    belowEvalThreshold: number;
}

export interface HHblitsHHInfoResult extends HHInfoResult {
    tm: number;
    coil: number;
    signal: string;
}

export interface HMMERHHInfoResult extends HHInfoResult {
    tm: number;
    coil: number;
    signal: string;
}

export interface HHompHHInfoResult extends HHInfoResult {
    probOMP: string;
}

export interface HHpredHHInfoResult extends HHInfoResult {
    tm: number;
    coil: number;
    signal: string;
    db: string;
    proteomes: number;
    modeller: boolean;
    msa_gen: string;
    qa3m_count: string;
    alignmentHash: string;
    query_neff: number;
}

export interface StructureFileResponse {
    data: any;
    filename?: string;
}

export interface ForwardingSubmission {
    forwardHitsMode: string;
    sequenceLengthMode: string;
    eval: number | string;
    selected: number[];
}
