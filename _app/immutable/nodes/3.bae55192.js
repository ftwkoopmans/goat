import{s as Le,f as a,a as l,l as Z,A as Se,g as i,d as t,c as r,h as q,m as ee,B as g,j as d,C as De,k as c,x as o,i as n,y as he,z as Ne}from"../chunks/scheduler.cd38125d.js";import{S as Ve,i as je}from"../chunks/index.a6a2bf8f.js";/* empty css                    */import{p as Fe}from"../chunks/stores.bacf73c7.js";function Be(p){let m,S,u,G,T,x,ge,te,I,se,ae,C,ie,oe,E,ne,le,O,re,B,H,pe="Documentation",R,D,fe=`<h3>Citation</h3> <div><p>The GOAT algorithm has not been published yet but a preprint is available, please cite it when
			using the early-access version of GOAT;
			<br/> <i>Koopmans, F. (2023). GOAT: efficient and robust identification of gene set enrichment.</i> <br/><a href="https://doi.org/10.1101/2023.12.10.570979" target="_blank" rel="nofollow">https://doi.org/10.1101/2023.12.10.570979</a></p></div>`,P,A,me=`<h3>How do I use this tool?</h3> <div><div>Overview of the GOAT online workflow</div> <div><ol><li class="svelte-pooqga">Input data: the gene list, this is the dataset you want to analyze. Typically a table of
					gene identifiers and their respective effect sizes that indicate an association with some
					experimental condition (e.g. summary statistics from an OMICS-based study). <br/>See
					below for the
					<a href="#genelist">Expected file format for your dataset</a> (also contains a link to an example
					dataset/gene list)</li> <li class="svelte-pooqga">Select the gene set database to perform enrichment testing on (e.g. the Gene Ontology
					database)</li> <li class="svelte-pooqga">Optionally, adjust the GOAT online settings</li> <li class="svelte-pooqga">run the analysis by clicking &quot;START&quot; (&quot;GOAT analysis&quot; section of the tool)</li> <li class="svelte-pooqga">Download results; an Excel table with all gene set statistics &amp; text file with a Methods
					description for your paper. <br/>See the
					<a href="#GOATresults">GOAT online output files</a> section below for a detailed description</li> <li class="svelte-pooqga">Use interactive data analysis tools to visualize/interpret results.</li></ol></div> <div>There is a <a href="#Glossary">Glossary</a> at the bottom of this page.</div></div> <br/> <b>Can I use the GOAT algorithm programmatically?</b> <div style="margin-top: 5px;">Yes! We also provide an R package; <a href="https://github.com/ftwkoopmans/goat" target="_blank" rel="nofollow">click here to go to the GitHub page</a></div>`,W,v,ve=`<h3>Expected file format for the input gene list:</h3> <div><ul><li class="svelte-pooqga">File format: either CSV, TSV, or Excel (.xlsx file, data on the first sheet). Note that for
				Excel, the old .xls format is not supported, only .xlsx files work with GOAT online.</li> <li class="svelte-pooqga">Required columns (column names must match exactly)</li> <ol><li class="svelte-pooqga"><b>gene</b>: Human Entrez (NCBI) gene identifiers (integer values)</li> <li class="svelte-pooqga"><b>symbol</b>: gene symbol (at least 2 characters)</li> <li class="svelte-pooqga"><b>effectsize</b>: effect size or log2 foldchange (numeric/decimal values)</li> <li class="svelte-pooqga"><b>pvalue</b>: gene p-values (numeric/decimal values). Use the un-adjusted p-values
					because after multiple comparisons adjustment, there typically are many more ties (e.g.
					p-values set to 1) among adjusted p-values which in turn causes a loss of information</li> <li class="svelte-pooqga"><b>signif</b>: was the adjusted p-value significant? This column should contain boolean
					values (true and false, or 0 and 1). Here one should use the adjusted p-values! While this
					information is NOT used by the GOAT algorithm to identify significant gene sets, flagging
					proteins that are significant in your gene list/dataset will yield useful information in
					the downstream interpretation of your data. For example, in the GOAT result tables, you can see
					for each gene set how many (and which) significant genes are present.
					<br/></li></ol> <li class="svelte-pooqga">Missing/empty values are not allowed, except for the &#39;gene&#39; column; rows where this column
				is empty will be ignored/skipped.</li> <li class="svelte-pooqga">Duplicate entries for genes, i.e. multiple rows with the same value in the &#39;gene&#39; column,
				are reduced to only 1 row; whichever has the lowest p-value (and if there is no p-value
				column, whichever row has the highest absolute effect size).</li> <li class="svelte-pooqga">While you can upload gene lists that only have either an &#39;effectsize&#39; or &#39;pvalue&#39; column, it
				is recommended to always include both if this information is available in your dataset
				because both sources of information can be used when sorting/ranking the gene list (e.g.
				when sorting by p-value, effect sizes can be used to break ties, especially for proteins with
				p-value=1).</li></ul></div> <div style="margin-bottom: 20px;"><b>importantly</b>, only Human Entrez (NCBI) gene IDs are supported for now (i.e. values in the
		&#39;gene&#39; column).
		<br/>
		We provide a gene ID mapping tool (available through the menu on top of this screen) to easily add
		Entrez gene IDs to your gene lists in case these only contain gene symbols.</div> <div style="margin-bottom: 10px;">Example gene list: <a href="data/example_genelist_Wingo_pmid32424284.xlsx">click here to download</a>
		the Wingo et al. 2020 dataset in Excel format (<a href="https://pubmed.ncbi.nlm.nih.gov/32424284" target="_blank" rel="nofollow">PMID:32424284</a>). You may use this as an example for preparing your gene list in a format compatible with
		GOAT.</div>`,K,b,be=`<h3>Using a custom gene set database</h3> <p>A gene set database prepared in the generic GMT file format can be easily imported into this tool.
		Any compatible GMT file that is stored on your computer can be used with GOAT online by using
		the &quot;upload GMT file&quot; button in the &quot;Gene sets&quot; section of the GOAT online tool.</p> <p>GOAT online only works when the input gene list and gene set database use Entrez Gene
		identifiers, so it is crucial to select the appropriate gene format when downloading a gene set
		database.</p> <p>The MSigDB collection website contains a large number of gene set databases that can be
		downloaded in GMT format and subsequently used in GOAT online;
		<a href="https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp" target="_blank" rel="nofollow">click here to visit their website</a>. To retrieve files compatible with the GOAT online tool, make sure to download via links that
		are labeled; <u>NCBI (Entrez) Gene IDs</u> (i.e. these are GMT files that contain collections of
		gene sets, with human gene identifiers in NCBI Entrez format).</p> <div>Example gene set databases that can be obtained from MSigDB (there are many more!):
		<ul><li class="svelte-pooqga">KEGG_MEDICUS: Canonical Pathways gene sets derived from the KEGG MEDICUS pathway database</li> <li class="svelte-pooqga">HPO: Gene sets derived from the Human Phenotype ontology</li> <li class="svelte-pooqga">CP: Canonical pathways. For example, the downloaded file would have a name similar to
				&quot;c2.cp.v2023.2.Hs.entrez.gmt&quot;</li> <li class="svelte-pooqga">C3: regulatory target gene sets. Gene sets representing potential targets of regulation by
				transcription factors or microRNAs</li></ul></div>`,U,y,ye=`<h3>GOAT online output files</h3> <p>After the GOAT analysis is done, you can download the results from the &quot;Result summary&quot; section of the
		GOAT online tool.</p> <p>The downloaded ZIP file contains a plain text file that details the exact settings used in GOAT
		online and includes a paragraph that can be used as Methods text in your manuscript (i.e. GOAT
		online version, citation, gene set database that was used, settings for multiple testing
		correction, etcetera).</p> <p>The Download also includes an Excel table that describes all tested gene sets and their
		GOAT-estimated p-values. These results can also be browsed in the online tool: click on &quot;TABLE&quot;
		in the &quot;Data analysis&quot; section.</p> <div>The following list details all columns that are provided in the Excel table. For data that is
		also shown in the &quot;TABLE&quot; section of the online tool, the respective column names are shown
		between brackets:
		<ul><li class="svelte-pooqga">source: for each gene set this describes the source/domain/classification that was used to
				categorize the gene set in the input gene set database. For example, in the GO database
				possible values for &quot;source&quot; are GO_CC, GO_BP, GO_MF describing the respective ontology
				domains Cellular Component, Biological Process, Molecular Function.</li> <li class="svelte-pooqga">id: the gene set identifier provided by the input gene set database</li> <li class="svelte-pooqga">name [gene set name]: the gene set &quot;name&quot; provided by the input gene set database</li> <li class="svelte-pooqga">ngenes_input: number of (unique) genes in the gene set as provided in the input gene set
				database</li> <li class="svelte-pooqga">ngenes [#genes]: number of genes that overlap between the input gene set and your gene list</li> <li class="svelte-pooqga">ngenes_signif [#signif]: number of genes that overlap between the input gene set and those
				genes that are flagged as &#39;signif&#39; in your gene list</li> <li class="svelte-pooqga">genes: gene symbols for the top 250 gene constituents (sorted by absolute effect size or
				p-value, whichever was used as &#39;score_type&#39;). The symbols are retrieved from your input gene
				list.</li> <li class="svelte-pooqga">genes_signif: gene symbols for the top 250 gene constituents that are flagged as &#39;signif&#39; in
				your input gene list (sorted by absolute effect size or p-value, whichever was used as
				&#39;score_type&#39;). The symbols are retrieved from your input gene list.</li> <li class="svelte-pooqga">zscore: a standardized z-score is computed from gene set p-values + effect size direction
				(up/down) if tested. Importantly, we here return standardized z-scores because the GOAT
				gene set score (mean of gene scores) is relative to the respective gene set-size-matched null
				distributions (a skewed normal)! In contrast, the standardized z-scores are comparable
				between gene sets (as are the p-values).
				<br/>
				Negative z-scores are returned for gene sets with score_type &#39;effectsize_down&#39;. i.e. these are
				gene sets that were explicitly tested for overrepresentation of negative effect size values.
				Note that &#39;effectsize_abs&#39; and &#39;pvalue&#39; score types are agnostic to up/down regulation so z-scores
				are always positive in that case.</li> <li class="svelte-pooqga">pvalue [p-value]: GOAT-estimated gene set p-value (as-is)</li> <li class="svelte-pooqga">pvalue_adjust [adj. p-value]: GOAT-estimated gene set p-value, after multiple testing
				corrections (depending on your chosen settings)</li> <li class="svelte-pooqga">signif [signif]: boolean value indicating whether the gene set was significant under your
				specified criteria (type of correction and cutoff value)</li> <li class="svelte-pooqga">score_type [score type]: the values shown here depend on the &quot;gene score type&quot; option that
				you selected in the &quot;GOAT analysis&quot; section of the web tool. If gene sets were tested by type
				&quot;effectsize&quot;, then GOAT will test for each gene set whether enrichment is strongest in
				either up- or down-regulation; if &quot;effectsize_up&quot;, testing the gene set for enrichment in
				positive gene effect sizes yields a lower (stronger) p-value than testing the gene set in
				the direction of negative gene effect sizes. Gene sets indicated as &quot;effectsize_down&quot; are the
				opposite; their gene constituents are more enriched for negative effect sizes. In case any
				other &quot;gene score type&quot; than &quot;effectsize&quot; was selected, this value will default to the
				configured &quot;gene score type&quot;. Note that in the online tool, &quot;effectsize&quot; is abbreviated to
				&quot;ES&quot;.</li></ul></div>`,Y,h,N,we="Description of the GOAT algorithm",ue,V,qe=`In brief, the Gene set Ordinal Association Test (GOAT) is a parameter-free permutation-based
		algorithm for gene set enrichment analysis. It is easy to use via the online web tool or R
		package, computationally efficient, and the resulting gene set p-values are well calibrated under
		the null hypothesis and invariant to gene set size. Application to various real-world proteomics
		and gene expression studies demonstrates that GOAT consistently identifies more significant Gene
		Ontology terms as compared to alternative methods.`,de,F,j,Ge,ce,L,Te=`<ol><li class="svelte-pooqga">Required input is a list of genes and their respective test statistics (p-value/effect
				size), and a gene set database obtained from GO or alternative resources.</li> <li class="svelte-pooqga">Test statistics from the gene list are transformed to gene scores by rank(-pvalue)^2 or
				rank(effect size)^2 depending on user input, i.e. smaller p-values translate to higher gene
				scores. The result is a skewed gene score distribution.</li> <li class="svelte-pooqga">For each gene set size N (number of genes), bootstrapping procedures generate a null
				distribution of gene set scores. This yields a skew-normal distribution for small gene sets
				and converges to a normal distribution for large gene sets.</li> <li class="svelte-pooqga">Gene set significance is determined for each gene set by comparing its score (mean of
				respective gene scores) against a null distribution of the same size (N).</li></ol>`,$,_,xe=`<h3>Frequently Asked Questions (FAQ)</h3> <p><i>What if I want to obtain the significant gene sets that also contain at least N &quot;significant
			genes&quot; from my input gene list/dataset?</i> <br/>The output table from this tool contains a column &quot;ngenes_signif&quot;, representing the number
		of genes that overlap between the gene set database and those genes that are flagged as &#39;signif&#39;
		in your gene list. You can use the information in this column for post-hoc filtering.</p> <p><i>How do we future-proof this tool? e.g. prevent it from growing stale/outdated and keep it
			online?</i> <br/>The website is implemented as a fully &quot;static&quot; HTML + Javascript website, so hosting the
		website is trivial; we use GitHub Pages to host it. As long as GitHub is online, so is this
		tool.
		<br/>To ensure access to recent Gene Ontology database data, we are setting up an automated
		workflow for importing the latest GO database biannually using GitHub Actions (i.e. does not
		require any manual updates/edits to this website). In the next major website update, you will
		be able to select from available GO database versions/snapshots. The current GO version
		available in GOAT online is 2024-01-01.
		<br/>Further, you can always import any gene set collection in GMT format into this web tool.</p> <p><i>Does the web tool yield the same gene set p-values as the R package?</i> <br/>Yes. We validated that the gene set p-values computed by GOAT online and the GOAT R
		package are the same across all OMICS-based datasets that are described in the GOAT manuscript.</p>`,J,k,Oe='<h3>Privacy</h3> <ul><li class="svelte-pooqga">all analyses are performed locally on your computer using client-side Javascript code</li> <li class="svelte-pooqga">your gene list and all analyses thereof remain private; your data does not leave your computer</li> <li class="svelte-pooqga">we do count the number of times this tool is used, anonymously, to gauge its popularity</li></ul>',Q,z,Ae=`<h3>Logo</h3> <p>the GOAT logo shown on this website is borrowed from the open-source Noto Emoji Font version
		14.0</p>`,X,w,_e=`<h3>Glossary</h3> <div><ul><li class="svelte-pooqga">GOAT: Gene set Ordinal Association Test. A parameter-free algorithm for gene set enrichment
				analysis of preranked gene lists.</li> <li class="svelte-pooqga">Gene list: A preranked gene list is here defined as a table of gene identifiers and their
				respective effect sizes and/or p-values that indicate association with some experimental
				condition (e.g. summary statistics from an OMICS-based study)</li> <li class="svelte-pooqga">Gene set: A gene set can be any set of genes of interest; it is typically defined as a set
				of genes that are known members of the same biological pathway, localized to the same
				(sub)cellular compartment, co-expressed under certain conditions or associated with some
				disorder as defined in a gene set database such as GO. In the GOAT R package (and online
				tool) one tests for enrichment of top-ranked genes in the input gene list against each gene
				set from some collection/database</li> <li class="svelte-pooqga">GO: the Gene Ontology database. <i>&quot;The Gene Ontology (GO) knowledgebase is the world’s largest source of information on the
					functions of genes. This knowledge is both human-readable and machine-readable, and is a
					foundation for computational analysis of large-scale molecular biology and genetics
					experiments in biomedical research.&quot;</i>
				Reference:
				<a href="https://www.geneontology.org" target="_blank" rel="nofollow">www.geneontology.org</a></li> <li class="svelte-pooqga">SynGO: the Synaptic Gene Ontology database. <i>&quot;An evidence-based, expert-curated resource for synapse function and gene enrichment
					studies&quot;</i>
				Reference:
				<a href="https://www.syngoportal.org" target="_blank" rel="nofollow">www.syngoportal.org</a></li> <li class="svelte-pooqga">KEGG: Kyoto Encyclopedia of Genes and Genomes. <i>&quot;KEGG is a database resource for understanding high-level functions and utilities of the
					biological system, such as the cell, the organism and the ecosystem, from molecular-level
					information, especially large-scale molecular datasets generated by genome sequencing and
					other high-throughput experimental technologies.&quot;</i>
				Reference:
				<a href="https://www.genome.jp/kegg/" target="_blank" rel="nofollow">www.genome.jp/kegg/</a>. Note that the easiest way to use KEGG pathways in GOAT online is to download respective
				GMT files from MSigDB, as detailed in the &quot;Using a custom gene set database&quot; section of this
				page.</li> <li class="svelte-pooqga">M&amp;M / Methods: The Materials and Methods section of a (scientific) manuscript, where one
				would typically detail exactly how GOAT was used to generate presented results.
				Note that results from GOAT include ready-made text that can be used for this.</li></ul></div>`;return{c(){m=a("meta"),S=l(),u=a("nav"),G=a("div"),T=a("a"),x=a("img"),te=l(),I=a("a"),se=Z("Home"),ae=l(),C=a("a"),ie=Z("GOAT online"),oe=l(),E=a("a"),ne=Z("gene ID mapping"),le=l(),O=a("a"),re=Z("Documentation"),B=l(),H=a("h1"),H.textContent=pe,R=l(),D=a("div"),D.innerHTML=fe,P=l(),A=a("div"),A.innerHTML=me,W=l(),v=a("div"),v.innerHTML=ve,K=l(),b=a("div"),b.innerHTML=be,U=l(),y=a("div"),y.innerHTML=ye,Y=l(),h=a("div"),N=a("h3"),N.textContent=we,ue=l(),V=a("p"),V.textContent=qe,de=l(),F=a("div"),j=a("img"),ce=l(),L=a("div"),L.innerHTML=Te,$=l(),_=a("div"),_.innerHTML=xe,J=l(),k=a("div"),k.innerHTML=Oe,Q=l(),z=a("div"),z.innerHTML=Ae,X=l(),w=a("div"),w.innerHTML=_e,this.h()},l(e){const s=Se("svelte-1nmuhnx",document.head);m=i(s,"META",{name:!0,content:!0}),s.forEach(t),S=r(e),u=i(e,"NAV",{});var f=q(u);G=i(f,"DIV",{style:!0});var ke=q(G);T=i(ke,"A",{href:!0});var ze=q(T);x=i(ze,"IMG",{src:!0,width:!0,height:!0,alt:!0}),ze.forEach(t),ke.forEach(t),te=r(f),I=i(f,"A",{href:!0,style:!0});var Me=q(I);se=ee(Me,"Home"),Me.forEach(t),ae=r(f),C=i(f,"A",{href:!0,style:!0});var Ie=q(C);ie=ee(Ie,"GOAT online"),Ie.forEach(t),oe=r(f),E=i(f,"A",{href:!0,style:!0});var Ce=q(E);ne=ee(Ce,"gene ID mapping"),Ce.forEach(t),le=r(f),O=i(f,"A",{href:!0,style:!0,class:!0});var Ee=q(O);re=ee(Ee,"Documentation"),Ee.forEach(t),f.forEach(t),B=r(e),H=i(e,"H1",{"data-svelte-h":!0}),g(H)!=="svelte-fjilyk"&&(H.textContent=pe),R=r(e),D=i(e,"DIV",{"data-svelte-h":!0}),g(D)!=="svelte-1g9y5mm"&&(D.innerHTML=fe),P=r(e),A=i(e,"DIV",{style:!0,"data-svelte-h":!0}),g(A)!=="svelte-uci5mp"&&(A.innerHTML=me),W=r(e),v=i(e,"DIV",{id:!0,style:!0,"data-svelte-h":!0}),g(v)!=="svelte-10hgk2c"&&(v.innerHTML=ve),K=r(e),b=i(e,"DIV",{id:!0,style:!0,"data-svelte-h":!0}),g(b)!=="svelte-krlbem"&&(b.innerHTML=be),U=r(e),y=i(e,"DIV",{id:!0,style:!0,"data-svelte-h":!0}),g(y)!=="svelte-k6u1pa"&&(y.innerHTML=ye),Y=r(e),h=i(e,"DIV",{style:!0});var M=q(h);N=i(M,"H3",{"data-svelte-h":!0}),g(N)!=="svelte-1xztreg"&&(N.textContent=we),ue=r(M),V=i(M,"P",{"data-svelte-h":!0}),g(V)!=="svelte-11qwqws"&&(V.textContent=qe),de=r(M),F=i(M,"DIV",{});var He=q(F);j=i(He,"IMG",{src:!0,alt:!0}),He.forEach(t),ce=r(M),L=i(M,"DIV",{style:!0,"data-svelte-h":!0}),g(L)!=="svelte-n2xzyy"&&(L.innerHTML=Te),M.forEach(t),$=r(e),_=i(e,"DIV",{style:!0,"data-svelte-h":!0}),g(_)!=="svelte-d07yls"&&(_.innerHTML=xe),J=r(e),k=i(e,"DIV",{style:!0,"data-svelte-h":!0}),g(k)!=="svelte-ro74ze"&&(k.innerHTML=Oe),Q=r(e),z=i(e,"DIV",{style:!0,"data-svelte-h":!0}),g(z)!=="svelte-curqo9"&&(z.innerHTML=Ae),X=r(e),w=i(e,"DIV",{id:!0,style:!0,"data-svelte-h":!0}),g(w)!=="svelte-s4pm5i"&&(w.innerHTML=_e),this.h()},h(){document.title="GOAT: Gene set Ordinal Association Test - documentation",d(m,"name","description"),d(m,"content","Gene set enrichment analysis for Gene Ontology (GO) or KEGG pathways using the GOAT algorithm web tool. Online data analysis for your preranked gene list from e.g. proteomics or bulk/scRNAseq gene expression studies"),De(x.src,ge=p[0]+"android-chrome-192x192.png")||d(x,"src",ge),d(x,"width","40"),d(x,"height","40"),d(x,"alt","GOAT"),d(T,"href",p[0]),c(G,"padding","4px"),c(G,"margin-left","20px"),d(I,"href",p[0]),c(I,"margin-left","5px"),d(C,"href",p[0]+"goat"),c(C,"margin-left","40px"),d(E,"href",p[0]+"genemap"),c(E,"margin-left","40px"),d(O,"href",p[0]+"docs"),c(O,"margin-left","40px"),d(O,"class","active"),c(A,"margin-top","50px"),d(v,"id","genelist"),c(v,"margin-top","50px"),d(b,"id","GMTgenesets"),c(b,"margin-top","50px"),d(y,"id","GOATresults"),c(y,"margin-top","50px"),De(j.src,Ge=p[0]+"goat_algorithm.png")||d(j,"src",Ge),d(j,"alt","GOAT algorithm"),c(L,"padding","0px 50px 0px 10px"),c(h,"margin-top","50px"),c(_,"margin-top","50px"),c(k,"margin-top","50px"),c(z,"margin-top","50px"),d(w,"id","Glossary"),c(w,"margin-top","50px")},m(e,s){o(document.head,m),n(e,S,s),n(e,u,s),o(u,G),o(G,T),o(T,x),o(u,te),o(u,I),o(I,se),o(u,ae),o(u,C),o(C,ie),o(u,oe),o(u,E),o(E,ne),o(u,le),o(u,O),o(O,re),n(e,B,s),n(e,H,s),n(e,R,s),n(e,D,s),n(e,P,s),n(e,A,s),n(e,W,s),n(e,v,s),n(e,K,s),n(e,b,s),n(e,U,s),n(e,y,s),n(e,Y,s),n(e,h,s),o(h,N),o(h,ue),o(h,V),o(h,de),o(h,F),o(F,j),o(h,ce),o(h,L),n(e,$,s),n(e,_,s),n(e,J,s),n(e,k,s),n(e,Q,s),n(e,z,s),n(e,X,s),n(e,w,s)},p:he,i:he,o:he,d(e){e&&(t(S),t(u),t(B),t(H),t(R),t(D),t(P),t(A),t(W),t(v),t(K),t(b),t(U),t(y),t(Y),t(h),t($),t(_),t(J),t(k),t(Q),t(z),t(X),t(w)),t(m)}}}function Re(p,m,S){let u;return Ne(p,Fe,T=>S(1,u=T)),[!!u&&u.url.hostname.endsWith("github.io")&&"/goat/"||"/"]}class Ye extends Ve{constructor(m){super(),je(this,m,Re,Be,Le,{})}}export{Ye as component};
