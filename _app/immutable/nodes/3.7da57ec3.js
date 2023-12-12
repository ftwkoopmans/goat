import{s as le,f as o,a as l,A as re,g as n,d as t,c as r,h as ne,B as u,j as k,D as O,k as p,x as d,i as s,y as S,z as ue}from"../chunks/scheduler.a30ae394.js";import{S as pe,i as he}from"../chunks/index.49b3fccf.js";/* empty css                    */import{p as de}from"../chunks/stores.b6b25a96.js";import{e as w}from"../chunks/singletons.3d14f100.js";function ce(x){let m,I,a,h,g,D=`<img src="${w+"/android-chrome-192x192.png"}" width="40" height="40" alt="GOAT"/>`,K,v,$="Home",R,f,Y="GOAT online",W,y,J="gene ID mapping",B,b,Q="Documentation",E,M,U="Documentation",L,H,Z=`<h3>Citation</h3> <div><p>The GOAT algorithm has not been published yet but a preprint is available, please cite it when
			using the early-access version of GOAT;
			<br/> <i>Koopmans, F. (2023). GOAT: efficient and robust identification of geneset enrichment.</i> <br/><a href="https://doi.org/10.1101/2023.12.10.570979" target="_blank" rel="nofollow">https://doi.org/10.1101/2023.12.10.570979</a></p> <p>A ready-made M&amp;M text describing your GOAT analysis is included with this tool&#39;s results
			(c.f. the log file contained in the output ZIP-file).</p></div>`,P,T,X='<h3>How do I use this tool?</h3> <div><div>brief overview of workflow</div> <div><ol><li class="svelte-pooqga">input data: genelist</li> <li class="svelte-pooqga">input data: genesets</li> <li class="svelte-pooqga">settings</li> <li class="svelte-pooqga">start</li> <li class="svelte-pooqga">view summary table + download link (contains Excel table and M&amp;M for your paper)</li> <li class="svelte-pooqga">use interactive data analysis tools to inspect results</li></ol></div></div> <br/> <b>Can I use the GOAT algorithm programmatically?</b> <div style="margin-top: 5px;">Yes! We also provide a R package; <a href="https://github.com/ftwkoopmans/goat" target="_blank" rel="nofollow">via this link</a></div>',z,_,ee=`<h3>Expected file format for the input genelist:</h3> <div><ul><li class="svelte-pooqga">File format: either CSV, TSV or Excel (.xlsx file, data on the first sheet). Note that for
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
		in a format compatible with GOAT.</div> <div>Example genelist: <a href="Klaassen_2016_IP_mass-spec_PMID26931375.xlsx">click here to download</a> the Klaassen et al. (2016, PMID:26931375) dataset in Excel format, an additional example.</div>`,N,G,te=`<h3>Using custom genesets</h3> <div>Genesets prepared in the generic GMT file format can be imported. For example, KEGG pathways in
		GMT format can be downloaded <a href="https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp" target="_blank" rel="nofollow">via this link</a>. Navigate to the &quot;KEGG_LEGACY subset of CP&quot; section and select &quot;NCBI (Entrez) Gene IDs&quot;. This
		should yield a file with a name similar to &quot;c2.cp.kegg.v2023.1.Hs.entrez.gmt&quot; (depending on the
		version). GOAT online only works when input genelists and genesets use Entrez Gene identifiers,
		so it is crucial to select the appropriate gene format when preparing/downloading genesets in
		GMT format.
		<br/>You may use the genesets in this file using the &quot;upload GMT file&quot; button in the &quot;Genesets&quot;
		section of the GOAT online tool.</div>`,V,q,ie=`<h3>Frequently Asked Questions (FAQ)</h3> <p><i>How do we future-proof this tool? e.g. prevent it from growing stale/outdated and keep it
			online ?</i> <br/>The website is implemented as a fully &quot;static&quot; HTML + Javascript website, so hosting the
		website is trivial; we use GitHub Pages to host it. As long as GitHub is online, so is this
		tool.
		<br/>To ensure access to recent Gene Ontology database data, we are currently setting up an
		automated workflow for importing the latest GO database N times per year using GitHub Actions
		(i.e. won&#39;t rely on us manually updating this website).
		<br/>Further, you can always import any geneset collection in GMT format into this webtool.</p> <p><i>Does the webtool yield different geneset p-values than the R package?</i> <br/>No. We validated that the geneset p-values computed by GOAT online and the GOAT R package
		are the same across all 7 &#39;omics datasets that are described in the GOAT manuscript.</p>`,j,A,ae='<h3>Privacy</h3> <ul><li class="svelte-pooqga">all analyses are performed locally on your computer using client-side Javascript code</li> <li class="svelte-pooqga">your genelist and all analyses thereof remain private, your data does not leave your computer</li> <li class="svelte-pooqga">we do keep a counter for the number of times this tool is used to gauge its popularity</li></ul>',F,C,se=`<h3>Logo</h3> <p>the GOAT logo shown on this website is borrowed from the open source Noto Emoji Font version
		14.0</p>`;return{c(){m=o("meta"),I=l(),a=o("nav"),h=o("div"),g=o("a"),g.innerHTML=D,K=l(),v=o("a"),v.textContent=$,R=l(),f=o("a"),f.textContent=Y,W=l(),y=o("a"),y.textContent=J,B=l(),b=o("a"),b.textContent=Q,E=l(),M=o("h1"),M.textContent=U,L=l(),H=o("div"),H.innerHTML=Z,P=l(),T=o("div"),T.innerHTML=X,z=l(),_=o("div"),_.innerHTML=ee,N=l(),G=o("div"),G.innerHTML=te,V=l(),q=o("div"),q.innerHTML=ie,j=l(),A=o("div"),A.innerHTML=ae,F=l(),C=o("div"),C.innerHTML=se,this.h()},l(e){const i=re("svelte-1m51m41",document.head);m=n(i,"META",{name:!0,content:!0}),i.forEach(t),I=r(e),a=n(e,"NAV",{});var c=ne(a);h=n(c,"DIV",{style:!0});var oe=ne(h);g=n(oe,"A",{href:!0,"data-svelte-h":!0}),u(g)!=="svelte-6nvu5r"&&(g.innerHTML=D),oe.forEach(t),K=r(c),v=n(c,"A",{href:!0,style:!0,"data-svelte-h":!0}),u(v)!=="svelte-xsa54j"&&(v.textContent=$),R=r(c),f=n(c,"A",{href:!0,style:!0,"data-svelte-h":!0}),u(f)!=="svelte-o349n8"&&(f.textContent=Y),W=r(c),y=n(c,"A",{href:!0,style:!0,"data-svelte-h":!0}),u(y)!=="svelte-sn5cge"&&(y.textContent=J),B=r(c),b=n(c,"A",{href:!0,style:!0,"data-svelte-h":!0}),u(b)!=="svelte-1dikf5g"&&(b.textContent=Q),c.forEach(t),E=r(e),M=n(e,"H1",{"data-svelte-h":!0}),u(M)!=="svelte-fjilyk"&&(M.textContent=U),L=r(e),H=n(e,"DIV",{"data-svelte-h":!0}),u(H)!=="svelte-1vomu16"&&(H.innerHTML=Z),P=r(e),T=n(e,"DIV",{style:!0,"data-svelte-h":!0}),u(T)!=="svelte-b1u3qb"&&(T.innerHTML=X),z=r(e),_=n(e,"DIV",{style:!0,"data-svelte-h":!0}),u(_)!=="svelte-14m738q"&&(_.innerHTML=ee),N=r(e),G=n(e,"DIV",{style:!0,"data-svelte-h":!0}),u(G)!=="svelte-ktgcbp"&&(G.innerHTML=te),V=r(e),q=n(e,"DIV",{style:!0,"data-svelte-h":!0}),u(q)!=="svelte-1wq36b7"&&(q.innerHTML=ie),j=r(e),A=n(e,"DIV",{style:!0,"data-svelte-h":!0}),u(A)!=="svelte-ixuru7"&&(A.innerHTML=ae),F=r(e),C=n(e,"DIV",{style:!0,"data-svelte-h":!0}),u(C)!=="svelte-8ovfie"&&(C.innerHTML=se),this.h()},h(){document.title="GOAT: Geneset Ordinal Association Test - documentation",k(m,"name","description"),k(m,"content","Geneset enrichment analysis for Geneset Ontology (GO) or KEGG pathways using the GOAT algorithm webtool. Online data analysis for your preranked genelist from e.g. proteomics or bulk/scRNAseq gene expression studies"),k(g,"href",w+"/"),O(g,"active",x[0]==="home"),p(h,"padding","4px"),p(h,"margin-left","20px"),k(v,"href",w+"/"),p(v,"margin-left","5px"),O(v,"active",x[0]==="home"),k(f,"href",w+"/goat"),p(f,"margin-left","40px"),O(f,"active",x[0]==="goat"),k(y,"href",w+"/genemap"),p(y,"margin-left","40px"),O(y,"active",x[0]==="genemap"),k(b,"href",w+"/docs"),p(b,"margin-left","40px"),O(b,"active",x[0]==="docs"),p(T,"margin-top","50px"),p(_,"margin-top","50px"),p(G,"margin-top","50px"),p(q,"margin-top","50px"),p(A,"margin-top","50px"),p(C,"margin-top","50px")},m(e,i){d(document.head,m),s(e,I,i),s(e,a,i),d(a,h),d(h,g),d(a,K),d(a,v),d(a,R),d(a,f),d(a,W),d(a,y),d(a,B),d(a,b),s(e,E,i),s(e,M,i),s(e,L,i),s(e,H,i),s(e,P,i),s(e,T,i),s(e,z,i),s(e,_,i),s(e,N,i),s(e,G,i),s(e,V,i),s(e,q,i),s(e,j,i),s(e,A,i),s(e,F,i),s(e,C,i)},p:S,i:S,o:S,d(e){e&&(t(I),t(a),t(E),t(M),t(L),t(H),t(P),t(T),t(z),t(_),t(N),t(G),t(V),t(q),t(j),t(A),t(F),t(C)),t(m)}}}function me(x,m,I){let a;ue(x,de,D=>I(1,a=D));const h=a&&a.url?a.url.pathname.replace(".html",""):"/";return[h===w+"/goat"&&"goat"||h===w+"/docs"&&"docs"||h===w+"/genemap"&&"genemap"||"home"]}class we extends pe{constructor(m){super(),he(this,m,me,ce,le,{})}}export{we as component};