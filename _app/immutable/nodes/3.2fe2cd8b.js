import{s as le,f as o,a as r,A as re,g as n,d as i,c as u,h as oe,B as p,j as I,C as y,k as h,x as d,i as a,y as ne,z as ue}from"../chunks/scheduler.7144be00.js";import{S as pe,i as he}from"../chunks/index.3da70c46.js";import{i as l}from"../chunks/util.14f89131.js";import{p as de}from"../chunks/stores.bbded9a9.js";function ce(c){let w,O,s,x,T,W='<img src="./android-chrome-192x192.png" width="40" height="40" alt="GOAT"/>',F,g,B="Home",K,m,Y="GOAT online",R,v,$="gene ID mapping",S,f,J="Documentation",D,H,Q="Documentation",E,k,U=`<h3>Citation</h3> <div><p>The GOAT algorithm has not been published yet but a preprint is available, please cite it when
			using the early-access version of GOAT;
			<br/> <i>Koopmans, F. (2023). GOAT: efficient and robust identification of geneset enrichment.</i> <br/><a href="https://doi.org/10.1101/2023.12.10.570979" target="_blank" rel="nofollow">https://doi.org/10.1101/2023.12.10.570979</a></p> <p>A ready-made M&amp;M text describing your GOAT analysis is included with this tool&#39;s results
			(c.f. the log file contained in the output ZIP-file).</p></div>`,L,_,Z='<h3>How do I use this tool?</h3> <div><div>brief overview of workflow</div> <div><ol><li class="svelte-pooqga">input data: genelist</li> <li class="svelte-pooqga">input data: genesets</li> <li class="svelte-pooqga">settings</li> <li class="svelte-pooqga">start</li> <li class="svelte-pooqga">view summary table + download link (contains Excel table and M&amp;M for your paper)</li> <li class="svelte-pooqga">use interactive data analysis tools to inspect results</li></ol></div></div> <br/> <b>Can I use the GOAT algorithm programmatically?</b> <div style="margin-top: 5px;">Yes! We also provide a R package; <a href="https://github.com/ftwkoopmans/goat" target="_blank" rel="nofollow">via this link</a></div>',z,G,X=`<h3>Expected file format for the input genelist:</h3> <div><ul><li class="svelte-pooqga">File format: either CSV, TSV or Excel (.xlsx file, data on the first sheet). Note that for
				Excel, the old .xls is not supported, only the newer .xlsx format.</li> <li class="svelte-pooqga">Required columns (column names must match exactly)</li> <ol><li class="svelte-pooqga"><b>gene</b>: Human Entrez (NCBI) gene IDs (integer values)</li> <li class="svelte-pooqga"><b>symbol</b>: gene symbol (at least 2 characters)</li> <li class="svelte-pooqga"><b>effectsize</b>: effectsize or log2 foldchange (numeric/decimal values)</li> <li class="svelte-pooqga"><b>pvalue</b>: gene p-values (numeric/decimal values). Preferably, use the un-adjusted
					p-values because this information is used for sorting and after adjustment, there
					typically are many more ties (e.g. p-values set to 1) which in turn causes a loss of
					information</li> <li class="svelte-pooqga"><b>signif</b>: was the adjusted p-value significant? This column should contain boolean
					values, i.e. (true and false, or 0 and 1). Here one should use the adjusted p-values!
					While this information is NOT used by the GOAT algorithm to identify significant genesets,
					flagging proteins that are significant in your genelist/dataset will yield useful
					information in downstream interpretation of your data. For example, in the GOAT result
					tables you can see for each geneset how many (and which) significant genes are present.
					<br/></li></ol> <li class="svelte-pooqga">Missing/empty values are not allowed, except for the &#39;gene&#39; column; rows where this column
				is empty will be ignored/skipped.</li> <li class="svelte-pooqga">Duplicate entries for genes, i.e. multiple rows with the same value in the &#39;gene&#39; column,
				are reduced to only 1 row; whichever has the lowest p-value (and if there is no p-value
				column, whichever row has the highest absolute effectsize).</li> <li class="svelte-pooqga">While you can upload genelists that only have either an &#39;effectsize&#39; or &#39;pvalue&#39; column, it
				is recommended to always include both if this information is available in your dataset
				because both sources of information can be used when sorting/ranking the genelist (e.g. when
				sorting by p-value, effectsizes can be used to break ties especially for proteins with
				p-value=1).</li></ul></div> <div style="margin-bottom: 20px;"><b>importantly</b>, only Human Entrez (NCBI) gene IDs are supported for now (i.e. values in the
		&#39;gene&#39; column).
		<br/>
		We provide a gene ID mapping tool (available through the menu on top of this screen) to easily add
		Entrez gene IDs to your genelists in case these only contain gene symbols.</div> <div style="margin-bottom: 10px;">Example genelist: <a href="Hondius_2021_mass-spec_PMID33492460.xlsx">click here to download</a> the
		Hondius et al. (2021, PMID:33492460) dataset in Excel format, an example for preparing your genelist
		in a format compatible with GOAT.</div> <div>Example genelist: <a href="Klaassen_2016_IP_mass-spec_PMID26931375.xlsx">click here to download</a> the Klaassen et al. (2016, PMID:26931375) dataset in Excel format, an additional example.</div>`,P,q,ee=`<h3>Using custom genesets</h3> <div>Genesets prepared in the generic GMT file format can be imported. For example, KEGG pathways in
		GMT format can be downloaded <a href="https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp" target="_blank" rel="nofollow">via this link</a>. Navigate to the &quot;KEGG_LEGACY subset of CP&quot; section and select &quot;NCBI (Entrez) Gene IDs&quot;. This
		should yield a file with a name similar to &quot;c2.cp.kegg.v2023.1.Hs.entrez.gmt&quot; (depending on the
		version). GOAT online only works when input genelists and genesets use Entrez Gene identifiers,
		so it is crucial to select the appropriate gene format when preparing/downloading genesets in
		GMT format.
		<br/>You may use the genesets in this file using the &quot;upload GMT file&quot; button in the &quot;Genesets&quot;
		section of the GOAT online tool.</div>`,N,A,te=`<h3>Frequently Asked Questions (FAQ)</h3> <p><i>How do we future-proof this tool? e.g. prevent it from growing stale/outdated and keep it
			online ?</i> <br/>The website is implemented as a fully &quot;static&quot; HTML + Javascript website, so hosting the
		website is trivial; we use GitHub Pages to host it. So as long as GitHub is online, so is this
		tool.
		<br/>To ensure access to recent Gene Ontology database data, we are currently setting up an
		automated workflow for importing the latest GO database N times per year using GitHub Actions
		(i.e. won&#39;t rely on us manually updating this website).
		<br/>Further, you can always import any geneset collection in GMT format into this webtool.</p> <p><i>Does the webtool yield different geneset p-values than the R package?</i> <br/>No. We validated that the geneset p-values computed by GOAT online and the GOAT R package
		are the same across all 7 &#39;omics datasets that are described in the GOAT manuscript.</p>`,V,C,ie='<h3>Privacy</h3> <ul><li class="svelte-pooqga">all analyses are performed locally on your computer using client-side Javascript code</li> <li class="svelte-pooqga">your genelist and all analyses thereof remain private, your data does not leave your computer</li> <li class="svelte-pooqga">we do keep a counter for the number of times this tool is used to gauge its popularity</li></ul>',j,M,ae=`<h3>Logo</h3> <p>the GOAT logo shown on this website is borrowed from the open source Noto Emoji Font version
		14.0</p>`;return{c(){w=o("meta"),O=r(),s=o("nav"),x=o("div"),T=o("a"),T.innerHTML=W,F=r(),g=o("a"),g.textContent=B,K=r(),m=o("a"),m.textContent=Y,R=r(),v=o("a"),v.textContent=$,S=r(),f=o("a"),f.textContent=J,D=r(),H=o("h1"),H.textContent=Q,E=r(),k=o("div"),k.innerHTML=U,L=r(),_=o("div"),_.innerHTML=Z,z=r(),G=o("div"),G.innerHTML=X,P=r(),q=o("div"),q.innerHTML=ee,N=r(),A=o("div"),A.innerHTML=te,V=r(),C=o("div"),C.innerHTML=ie,j=r(),M=o("div"),M.innerHTML=ae,this.h()},l(e){const t=re("svelte-1m51m41",document.head);w=n(t,"META",{name:!0,content:!0}),t.forEach(i),O=u(e),s=n(e,"NAV",{});var b=oe(s);x=n(b,"DIV",{style:!0});var se=oe(x);T=n(se,"A",{href:!0,"data-svelte-h":!0}),p(T)!=="svelte-1kp1n6l"&&(T.innerHTML=W),se.forEach(i),F=u(b),g=n(b,"A",{href:!0,style:!0,"data-svelte-h":!0}),p(g)!=="svelte-9ooqgz"&&(g.textContent=B),K=u(b),m=n(b,"A",{href:!0,style:!0,"data-svelte-h":!0}),p(m)!=="svelte-10g8z77"&&(m.textContent=Y),R=u(b),v=n(b,"A",{href:!0,style:!0,"data-svelte-h":!0}),p(v)!=="svelte-13lvetu"&&(v.textContent=$),S=u(b),f=n(b,"A",{href:!0,style:!0,"data-svelte-h":!0}),p(f)!=="svelte-1106bxz"&&(f.textContent=J),b.forEach(i),D=u(e),H=n(e,"H1",{"data-svelte-h":!0}),p(H)!=="svelte-fjilyk"&&(H.textContent=Q),E=u(e),k=n(e,"DIV",{"data-svelte-h":!0}),p(k)!=="svelte-1vomu16"&&(k.innerHTML=U),L=u(e),_=n(e,"DIV",{style:!0,"data-svelte-h":!0}),p(_)!=="svelte-b1u3qb"&&(_.innerHTML=Z),z=u(e),G=n(e,"DIV",{style:!0,"data-svelte-h":!0}),p(G)!=="svelte-14m738q"&&(G.innerHTML=X),P=u(e),q=n(e,"DIV",{style:!0,"data-svelte-h":!0}),p(q)!=="svelte-ktgcbp"&&(q.innerHTML=ee),N=u(e),A=n(e,"DIV",{style:!0,"data-svelte-h":!0}),p(A)!=="svelte-o44mj1"&&(A.innerHTML=te),V=u(e),C=n(e,"DIV",{style:!0,"data-svelte-h":!0}),p(C)!=="svelte-ixuru7"&&(C.innerHTML=ie),j=u(e),M=n(e,"DIV",{style:!0,"data-svelte-h":!0}),p(M)!=="svelte-8ovfie"&&(M.innerHTML=ae),this.h()},h(){document.title="GOAT: Geneset Ordinal Association Test - documentation",I(w,"name","description"),I(w,"content","Geneset enrichment analysis for Geneset Ontology (GO) or KEGG pathways using the GOAT algorithm webtool. Online data analysis for your preranked genelist from e.g. proteomics or bulk/scRNAseq gene expression studies"),I(T,"href","/"),y(T,"active",l(c[0],"github.io/goat")||l(c[0],"/index")),h(x,"padding","4px"),h(x,"margin-left","20px"),I(g,"href","/"),h(g,"margin-left","5px"),y(g,"active",l(c[0],"github.io/goat")||l(c[0],"/index")),I(m,"href","/goat"),h(m,"margin-left","40px"),y(m,"active",l(c[0],"/goat")),I(v,"href","/genemap"),h(v,"margin-left","40px"),y(v,"active",l(c[0],"/genemap")),I(f,"href","/docs"),h(f,"margin-left","40px"),y(f,"active",l(c[0],"/docs")),h(_,"margin-top","50px"),h(G,"margin-top","50px"),h(q,"margin-top","50px"),h(A,"margin-top","50px"),h(C,"margin-top","50px"),h(M,"margin-top","50px")},m(e,t){d(document.head,w),a(e,O,t),a(e,s,t),d(s,x),d(x,T),d(s,F),d(s,g),d(s,K),d(s,m),d(s,R),d(s,v),d(s,S),d(s,f),a(e,D,t),a(e,H,t),a(e,E,t),a(e,k,t),a(e,L,t),a(e,_,t),a(e,z,t),a(e,G,t),a(e,P,t),a(e,q,t),a(e,N,t),a(e,A,t),a(e,V,t),a(e,C,t),a(e,j,t),a(e,M,t)},p(e,[t]){t&1&&y(T,"active",l(e[0],"github.io/goat")||l(e[0],"/index")),t&1&&y(g,"active",l(e[0],"github.io/goat")||l(e[0],"/index")),t&1&&y(m,"active",l(e[0],"/goat")),t&1&&y(v,"active",l(e[0],"/genemap")),t&1&&y(f,"active",l(e[0],"/docs"))},i:ne,o:ne,d(e){e&&(i(O),i(s),i(D),i(H),i(E),i(k),i(L),i(_),i(z),i(G),i(P),i(q),i(N),i(A),i(V),i(C),i(j),i(M)),i(w)}}}function ge(c,w,O){let s;return ue(c,de,x=>O(0,s=x)),[s]}class ye extends pe{constructor(w){super(),he(this,w,ge,ce,le,{})}}export{ye as component};