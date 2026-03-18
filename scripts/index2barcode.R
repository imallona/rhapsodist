v1_linker1 = 'ACTGGCCTGCGA'
v1_linker2 = 'GGTAGCGGTGACA'

Enh_linker1 = 'GTGA'
Enh_linker2 = 'GACA'

Enh_5p_primer = "ACAGGAAACTCATGGTGCGT"

Enh_5p_linker1 = "AATG"
Enh_5p_linker2 = "CCAC"

Enh_inserts = c("", "A", "GT", "TCA")

Tso_capture_seq_Enh_EnhV2 = "TATGCGTAGTAGGTATG"
Tso_capture_seq_EnhV3 = "GTGGAGTCGTGATTATA"

A96_cell_key1 = c("GTCGCTATA","CTTGTACTA","CTTCACATA","ACACGCCGG","CGGTCCAGG","AATCGAATG","CCTAGTATA","ATTGGCTAA","AAGACATGC","AAGGCGATC",
                 "GTGTCCTTA","GGATTAGGA","ATGGATCCA","ACATAAGCG","AACTGTATT","ACCTTGCGG","CAGGTGTAG","AGGAGATTA","GCGATTACA","ACCGGATAG",
                 "CCACTTGGA","AGAGAAGTT","TAAGTTCGA","ACGGATATT","TGGCTCAGA","GAATCTGTA","ACCAAGGAC","AGTATCTGT","CACACACTA","ATTAAGTGC",
                 "AAGTAACCC","AAATCCTGT","CACATTGCA","GCACTGTCA","ATACTTAGG","GCAATCCGA","ACGCAATCA","GAGTATTAG","GACGGATTA","CAGCTGACA",
                 "CAACATATT","AACTTCTCC","CTATGAAAT","ATTATTACC","TACCGAGCA","TCTCTTCAA","TAAGCGTTA","GCCTTACAA","AGCACACAG","ACAGTTCCG",
                 "AGTAAAGCC","CAGTTTCAC","CGTTACTAA","TTGTTCCAA","AGAAGCACT","CAGCAAGAT","CAAACCGCC","CTAACTCGC","AATATTGGG","AGAACTTCC",
                 "CAAAGGCAC","AAGCTCAAC","TCCAGTCGA","AGCCATCAC","AACGAGAAG","CTACAGAAC","AGAGCTATG","GAGGATGGA","TGTACCTTA","ACACACAAA",
                 "TCAGGAGGA","GAGGTGCTA","ACCCTGACC","ACAAGGATC","ATCCCGGAG","TATGTGGCA","GCTGCCAAT","ATCAGAGCT","TCGAAGTGA","ATAGACGAG",
                 "AGCCCAATC","CAGAATCGT","ATCTCCACA","ACGAAAGGT","TAGCTTGTA","ACACGAGAT","AACCGCCTC","ATTTAGATG","CAAGCAAGC","CAAAGTGTG",
                 "GGCAAGCAA","GAGCCAATA","ATGTAATGG","CCTGAGCAA","GAGTACATT","TGCGATCTA")

A96_cell_key2 = c("TACAGGATA","CACCAGGTA","TGTGAAGAA","GATTCATCA","CACCCAAAG","CACAAAGGC","GTGTGTCGA","CTAGGTCCT","ACAGTGGTA","TCGTTAGCA",
                 "AGCGACACC","AAGCTACTT","TGTTCTCCA","ACGCGAAGC","CAGAAATCG","ACCAAAATG","AGTGTTGTC","TAGGGATAC","AGGGCTGGT","TCATCCTAA",
                 "AATCCTGAA","ATCCTAGGA","ACGACCACC","TTCCATTGA","TAGTCTTGA","ACTGTTAGA","ATTCATCGT","ACTTCGAGC","TTGCGTACA","CAGTGCCCG",
                 "GACACTTAA","AGGAGGCGC","GCCTGTTCA","GTACATCTA","AATCAGTTT","ACGATGAAT","TGACAGACA","ATTAGGCAT","GGAGTCTAA","TAGAACACA",
                 "AAATAAATA","CCGACAAGA","CACCTACCC","AAGAGTAGA","TCATTGAGA","GACCTTAGA","CAAGACCTA","GGAATGATA","AAACGTACC","ACTATCCTC",
                 "CCGTATCTA","ACACATGTC","TTGGTATGA","GTGCAGTAA","AGGATTCAA","AGAATGGAG","CTCTCTCAA","GCTAACTCA","ATCAACCGA","ATGAGTTAC",
                 "ACTTGATGA","ACTTTAACT","TTGGAGGTA","GCCAATGTA","ATCCAACCG","GATGAACTG","CCATGCACA","TAGTGACTA","AAACTGCGC","ATTACCAAG",
                 "CACTCGAGA","AACTCATTG","CTTGCTTCA","ACCTGAGTC","AGGTTCGCT","AAGGACTAT","CGTTCGGTA","AGATAGTTC","CAATTGATC","GCATGGCTA",
                 "ACCAGGTGT","AGCTGCCGT","TATAGCCCT","AGAGGACCA","ACAATATGG","CAGCACTTC","CACTTATGT","AGTGAAAGG","AACCCTCGG","AGGCAGCTA",
                 "AACCAAAGT","GAGTGCGAA","CGCTAAGCA","AATTATAAC","TACTAGTCA","CAACAACGG"
)

A96_cell_key3 = c("AAGCCTTCT","ATCATTCTG","CACAAGTAT","ACACCTTAG","GAACGACAA","AGTCTGTAC","AAATTACAG","GGCTACAGA","AATGTATCG","CAAGTAGAA",
                 "GATCTCTTA","AACAACGCG","GGTGAGTTA","CAGGGAGGG","TCCGTCTTA","TGCATAGTA","ACTTACGAT","TGTATGCGA","GCTCCTTGA","GGCACAACA",
                 "CTCAAGACA","ACGCTGTTG","ATATTGTAA","AAGTTTACG","CAGCCTGGC","CTATTAGCC","CAAACGTGG","AAAGTCATT","GTCTTGGCA","GATCAGCGA",
                 "ACATTCGGC","AGTAATTAG","TGAAGCCAA","TCTACGACA","CATAACGTT","ATGGGACTC","GATAGAGGA","CTACATGCG","CAACGATCT","GTTAGCCTA",
                 "AGTTGCATC","AAGGGAACT","ACTACATAT","CTAAGCTTC","ACGAACCAG","TACTTCGGA","AACATCCAT","AGCCTGGTT","CAAGTTTCC","CAGGCATTT",
                 "ACGTGGGAG","TCTCACGGA","GCAACATTA","ATGGTCCGT","CTATCATGA","CAATACAAG","AAAGAGGCC","GTAGAAGCA","GCTATGGAA","ACTCCAGGG",
                 "ACAAGTGCA","GATGGTCCA","TCCTCAATA","AATAAACAA","CTGTACGGA","CTAGATAGA","AGCTATGTG","AAATGGAGG","AGCCGCAAG","ACAGTAAAC",
                 "AACGTGTGA","ACTGAATTC","AAGGGTCAG","TGTCTATCA","TCAGATTCA","CACGATCCG","AACAGAAAC","CATGAATGA","CGTACTACG","TTCAGCTCA",
                 "AAGGCCGCA","GGTTGGACA","CGTCTAGGT","AATTCGGCG","CAACCTCCA","CAATAGGGT","ACAGGCTCC","ACAACTAGT","AGTTGTTCT","AATTACCGG",
                 "ACAAACTTT","TCTCGGTTA","ACTAGACCG","ACTCATACG","ATCGAGTCT","CATAGGTCA")


index_to_sequence <- function(index, bead_version) {
    
    zerobased <- as.integer(index) - 1
    
    cl1 <- (as.integer(zerobased / 384 / 384) %% 384) + 1
    cl2 <- (as.integer(zerobased / 384) %% 384) + 1
    cl3 <- (zerobased %% 384) + 1
    
    if (bead_version == "v1") {
        cls1_sequence <- A96_cell_key1[cl1]
        cls2_sequence <- A96_cell_key2[cl2]
        cls3_sequence <- A96_cell_key3[cl3]
        
        return(paste0(cls1_sequence, v1_linker1, cls2_sequence, 
                      v1_linker2, cls3_sequence))
        
    } else if (bead_version == "Enh") {
        
        diversityInsert <- ""
        
        if (cl1 >= 1 && cl1 <= 24) {
            diversityInsert <- ""
        } else if (cl1 >= 25 && cl1 <= 48) {
            diversityInsert <- "A"
        } else if (cl1 >= 49 && cl1 <= 72) {
            diversityInsert <- "GT"
        } else { # cl1 >= 73 && cl1 <= 96
            diversityInsert <- "TCA"
        }
        
        cls1_sequence <- A96_cell_key1[cl1]
        cls2_sequence <- A96_cell_key2[cl2]
        cls3_sequence <- A96_cell_key3[cl3]
        
        #return(paste0(diversityInsert, cls1_sequence, Enh_linker1, cls2_sequence, Enh_linker2, cls3_sequence))
        return(paste0(cls1_sequence, cls2_sequence, cls3_sequence))
        
        
    } else if (bead_version == "EnhV2") {

        cls1_sequence <- B384_cell_key1[cl1]
        cls2_sequence <- B384_cell_key2[cl2]
        cls3_sequence <- B384_cell_key3[cl3]

        return(paste0(cls1_sequence, cls2_sequence, cls3_sequence))
    }
}


