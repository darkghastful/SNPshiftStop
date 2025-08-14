SNPshiftStop: README
================
Bailey Quinn
2025-08-14

# **SNPshiftSearch**

------------------------------------------------------------------------

## **Installation**

<!-- From CRAN (when published) -->
<!-- install.packages("") -->

``` r
# install.packages("devtools")      
install.packages("darkghastful/SNPshiftStop")
```

------------------------------------------------------------------------

## **Quick Start**

### **Load the package**

``` r
library(SNPshiftStop)
```

### **Usage example**

``` r
CFTR.F405L <- SNPshiftStop(gene.symbol="CFTR", taxon="human", mutation.from.start=1213, ref.bp="T", ref.aa="F", alt.aa="L", aa.mut.pos=405)
```

------------------------------------------------------------------------

## **Functions**

### `SNPsearch(rsid)`

Uses an rsid to query NCBI and identifies all associated SNPs.

- **Arguments**
  - `rsid` — Rsid associated with SNP.
- **Returns**
  - A frame of SNPs.

### **Example**

#### `SNPsearch("rs333")`

    #>   gene  sequence_ontology        ref.dna                               ref.bp
    #> 1 CCR5               <NA>    NM_000579.4 TACAGTCAGTATCAATTCTGGAAGAATTTCCAGACA
    #> 2 CCR5               <NA> NM_001100168.2 TACAGTCAGTATCAATTCTGGAAGAATTTCCAGACA
    #> 3 CCR5               <NA> NM_001394783.1 TACAGTCAGTATCAATTCTGGAAGAATTTCCAGACA
    #> 4 CCR5 frameshift_variant    NM_000579.4    AGTCAGTATCAATTCTGGAAGAATTTCCAGACA
    #> 5 CCR5 frameshift_variant NM_001100168.2    AGTCAGTATCAATTCTGGAAGAATTTCCAGACA
    #> 6 CCR5 frameshift_variant NM_001394783.1    AGTCAGTATCAATTCTGGAAGAATTTCCAGACA
    #>                                 alt.bp        ref.pro
    #> 1 TACAGTCAGTATCAATTCTGGAAGAATTTCCAGACA    NP_000570.1
    #> 2 TACAGTCAGTATCAATTCTGGAAGAATTTCCAGACA NP_001093638.1
    #> 3 TACAGTCAGTATCAATTCTGGAAGAATTTCCAGACA NP_001381712.1
    #> 4                                    A    NP_000570.1
    #> 5                                    A NP_001093638.1
    #> 6                                    A NP_001381712.1
    #>                          hgvs
    #> 1      NM_000579.4:c.551_585=
    #> 2   NM_001100168.2:c.551_585=
    #> 3   NM_001394783.1:c.551_585=
    #> 4    NM_000579.4:c.554_585del
    #> 5 NM_001100168.2:c.554_585del
    #> 6 NM_001394783.1:c.554_585del

### `SNPshiftStop(gene.symbol, taxon="human", mutation.from.start=NA, ref.bp="-", alt.bp="-", ref.aa=NA, alt.aa=NA, aa.mut.pos=NA, start.loss=FALSE, ncbi.api.key=NULL)`

Provides information on an SNP.

- **Arguments**
  - `gene.symbol` - Gene symbol with frameshift mutation
  - `taxon` - provide taxon character string or number; found using NCBI
    (default human (9606))
  - `mutation.from.start` - Location of mutation relative to start or
  - `ref.bp` - reference base pair (default is ‘“-”’)
  - `alt.bp` - alternate base pair (default is ‘“-”’)
  - `ref.aa` - reference amino acid (default is NA)
  - `alt.aa` - alternate amino acid (default is NA)
  - `aa.mut.pos` - amino acid mutation position (default is NA)
  - `start.loss` - does this mutation contain a start.loss (default is
    TRUE)
- **Returns**  
  A summary of the results will be printed including the reference amino
  acid sequence and alternate amino acid sequence.
  - `gene.symbol` - Gene symbol.
  - `aa.seq` - Amino acid sequence of the protein.
  - `mutated.aa.seq` - Mutated amino acid sequence.
  - `shortened.of.total.aa` - Length of alternate amino acid
    sequence/length of origional.
  - `dna.seq` - DNA sequence encoding the protein.
  - `mutated.dna.seq` - Alternate DNA sequence.

### **Example**

#### `SNPshiftStop(gene.symbol="CFTR", taxon="human", mutation.from.start=1213, ref.bp="T", ref.aa="F", alt.aa="L", aa.mut.pos=405)`

    #> 
    #> [1mReference:[22m
    #>  MQRSPLEKASVVSKLFFSWTRPILRKGYRQRLELSDIYQIPSVDSADNLSEKLEREWDRELASKKNPKLINALRRCFFWRFMFYGIFLYLGEVTKAVQPLLLGRIIASYDPDNKEERSIAIYLGIGLCLLFIVRTLLLHPAIFGLHHIGMQMRIAMFSLIYKKTLKLSSRVLDKISIGQLVSLLSNNLNKFDEGLALAHFVWIAPLQVALLMGLIWELLQASAFCGLGFLIVLALFQAGLGRMMMKYRDQRAGKISERLVITSEMIENIQSVKAYCWEEAMEKMIENLRQTELKLTRKAAYVRYFNSSAFFFSGFFVVFLSVLPYALIKGIILRKIFTTISFCIVLRMAVTRQFPWAVQTWYDSLGAINKIQDFLQKQEYKTLEYNLTTTEVVMENVTAFWEEGFGELFEKAKQNNNNRKTSNGDDSLFFSNFSLLGTPVLKDINFKIERGQLLAVAGSTGAGKTSLLMVIMGELEPSEGKIKHSGRISFCSQFSWIMPGTIKENIIFGVSYDEYRYRSVIKACQLEEDISKFAEKDNIVLGEGGITLSGGQRARISLARAVYKDADLYLLDSPFGYLDVLTEKEIFESCVCKLMANKTRILVTSKMEHLKKADKILILHEGSSYFYGTFSELQNLQPDFSSKLMGCDSFDQFSAERRNSILTETLHRFSLEGDAPVSWTETKKQSFKQTGEFGEKRKNSILNPINSIRKFSIVQKTPLQMNGIEEDSDEPLERRLSLVPDSEQGEAILPRISVISTGPTLQARRRQSVLNLMTHSVNQGQNIHRKTTASTRKVSLAPQANLTELDIYSRRLSQETGLEISEEINEEDLKECFFDDMESIPAVTTWNTYLRYITVHKSLIFVLIWCLVIFLAEVAASLVVLWLLGNTPLQDKGNSTHSRNNSYAVIITSTSSYYVFYIYVGVADTLLAMGFFRGLPLVHTLITVSKILHHKMLHSVLQAPMSTLNTLKAGGILNRFSKDIAILDDLLPLTIFDFIQLLLIVIGAIAVVAVLQPYIFVATVPVIVAFIMLRAYFLQTSQQLKQLESEGRSPIFTHLVTSLKGLWTLRAFGRQPYFETLFHKALNLHTANWFLYLSTLRWFQMRIEMIFVIFFIAVTFISILTTGEGEGRVGIILTLAMNIMSTLQWAVNSSIDVDSLMRSVSRVFKFIDMPTEGKPTKSTKPYKNGQLSKVMIIENSHVKKDDIWPSGGQMTVKDLTAKYTEGGNAILENISFSISPGQRVGLLGRTGSGKSTLLSAFLRLLNTEGEIQIDGVSWDSITLQQWRKAFGVIPQKVFIFSGTFRKNLDPYEQWSDQEIWKVADEVGLRSVIEQFPGKLDFVLVDGGCVLSHGHKQLMCLARSVLSKAKILLLDEPSAHLDPVTYQIIRRTLKQAFADCTVILCEHRIEAMLECQQFLVIEENKVRQYDSIQKLLNERSLFRQAISPSDRVKLFPHRNSSKCKSKPQIAALKEETEEEVQDTRL 
    #> [1mAlternate:[22m
    #>  MQRSPLEKASVVSKLFFSWTRPILRKGYRQRLELSDIYQIPSVDSADNLSEKLEREWDRELASKKNPKLINALRRCFFWRFMFYGIFLYLGEVTKAVQPLLLGRIIASYDPDNKEERSIAIYLGIGLCLLFIVRTLLLHPAIFGLHHIGMQMRIAMFSLIYKKTLKLSSRVLDKISIGQLVSLLSNNLNKFDEGLALAHFVWIAPLQVALLMGLIWELLQASAFCGLGFLIVLALFQAGLGRMMMKYRDQRAGKISERLVITSEMIENIQSVKAYCWEEAMEKMIENLRQTELKLTRKAAYVRYFNSSAFFFSGFFVVFLSVLPYALIKGIILRKIFTTISFCIVLRMAVTRQFPWAVQTWYDSLGAINKIQDFLQKQEYKTLEYNLTTTEVVMENVTAFWEEGLGNYLRKQNKTITIEKLLMVMTASSSVISHFLVLLS*KILISR*KEDSCWRLLDPLEQARLHF*W*LWENWSLQRVKLSTVEEFHSVLSFPGLCLAPLKKISSLVFPMMNIDTEASSKHAN*KRTSPSLQRKTI*FLEKVESH*VEVNEQEFL*QEQYTKMLICIY*TLLLDT*MF*QKKKYLKAVSVN*WLTKLGFWSLLKWNI*RKLTKY*FCMKVAAIFMGHFQNSKIYSQTLAQNSWDVILSTNLVQKEEIQS*LRPYTVSH*KEMLLSPGQKQKNNLLNRLESLGKKGRILFSIQSTLYENFPLCKRLPYK*MASKRILMSL*REGCP*YQILSRERRYCLASA*SALAPRFRHEGGSLS*T**HTQLTKVRTFTERQQHPHEKCHWPLRQT*LNWIYIQEGYLKKLAWK*VKKLTKKT*RSAFLMIWRAYQQ*LHGTHTFDILLSTRA*FLC*FGA**FFWQRWLLLWLCCGSLETLLFKTKGIVLIVEITAMQ*LSPAPVRIMCFTFTWE*PTLCLLWDSSEVYHWCIL*SQCRKFYTTKCYILFFKHLCQPSTR*KQVGFLIDSPKI*QFWMTFCLLPYLTSSSCY*L*LEL*QLSQFYNPTSLLQQCQ**WLLLC*EHISSKPHSNSNNWNLKAGVQFSLILLQA*KDYGHFVPSDGSLTLKLCSTKL*IYILPTGSCTCQHCAGSK*E*K*FLSSSSLLLPSFPF*QQEKEKEELVLS*L*P*IS*VHCSGL*TPA*MWIA*CDL*AESLSSLTCQQKVNLPSQPNHTRMANSRKL*LLRIHT*RKMTSGPQGAK*LSKISQQNTQKVEMPY*RTFPSQ*VLARGWASWEELDQGRVLCYQLF*DY*TLKEKSRSMVCLGIQ*LCNSGGKPLE*YHRKYLFFLEHLEKTWIPMNSGVIKKYGKLQMRLGSDL**NSFLGSLTLSLWMGAVS*AMATSS*CAWLDLFSVRRRSCCLMNPVLIWIQ*HTK*LEEL*NKHLLIAQ*FSVNTG*KQCWNANNFWS*KRTKCGSTIPSRNC*TRGASSGKPSAPPTG*SSFPTGTQASASLSPRLLL*KRRQKKRCKIQGF 
    #> 
    #> [1mIn the gene CFTR there was an T deletion at position 1213. The frameshift began at 405 with a change of F to L. This introduced a premature stop at position 441/1480 (shortening the protein by 1039 amino acids).[22m

<!-- snp.query -->
<!-- API documentation[https://clinicaltables.nlm.nih.gov/apidoc/snps/v3/doc.html]: -->
<!-- Output for an API query is an array of the following elements: -->
<!-- - The total number of results on the server, which can be more than the number of results returned. This reported total number of results may also be significantly less than the actual number of results and is limited to 10,000, which may significantly improve the service response time. -->
<!-- - An array of codes for the returned items. (This is the field specified with the cf query parameter above.) -->
<!-- - A hash of the "extra" data requested via the "ef" query parameter above. The keys on the hash are the fields (or their requested aliases) named in the "ef" parameter, and the value for a field is an array of that field's values in the same order as the returned codes. -->
<!-- - An array, with one element for each returned code, where each element is an array of the display strings specified with the "df" query parameter. -->
<!-- - An array, with one element for each returned code, where each element is the "code system" for the returned code. Note that only code-system aware APIs will return this array. -->

------------------------------------------------------------------------

## **Dependencies**

- **bqutils**  
- **Biostrings**
- **stringr**  
- **magrittr**  
- **methods**

------------------------------------------------------------------------

## **License**

GPL-3.0 license © Bailey Quinn
