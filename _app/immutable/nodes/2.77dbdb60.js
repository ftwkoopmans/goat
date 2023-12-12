import{s as Xe,f as t,a as r,l as ee,A as Je,g as s,d as n,c as l,h,m as te,B as x,j as a,C as se,D as Q,k as d,x as e,i as X,y as Ae,z as Qe}from"../chunks/scheduler.a30ae394.js";import{S as We,i as Ye}from"../chunks/index.49b3fccf.js";/* empty css                    */import{p as Ze}from"../chunks/stores.b6b25a96.js";import{e as Ue}from"../chunks/singletons.3d14f100.js";function et(o){let b,B,i,p,f,y,W,ae,A,ne,ie,O,re,le,T,oe,de,C,ce,Y,q,Oe=`<h1 style="" class="svelte-1vmxgqy">GOAT online geneset enrichment analysis</h1> <div><h3 class="svelte-1vmxgqy">Identify enriched genesets in a preranked genelist generated by e.g. proteomics or gene
			expression studies using the GOAT algorithm.</h3></div>`,Z,g,M,Te=`<h2 class="svelte-1vmxgqy">The GOAT algorithm for geneset enrichment analysis</h2> <div><p>The GOAT algorithm has not been published yet but a preprint is available, please cite it
				when using the early-access version of GOAT;
				<br/> <i>Koopmans, F. (2023). GOAT: efficient and robust identification of geneset enrichment.</i> <br/> <a href="https://doi.org/10.1101/2023.12.10.570979" target="_blank" rel="nofollow">https://doi.org/10.1101/2023.12.10.570979</a></p></div> <div style="margin-top: 25px;">GOAT features:
			<ul class="svelte-1vmxgqy"><li class="svelte-1vmxgqy">Sensitive; more significant genesets as compared to other methods</li> <li class="svelte-1vmxgqy">Accurate; geneset p-values are accurate under the null hypothesis</li> <li class="svelte-1vmxgqy">Fast; completes in seconds</li> <li class="svelte-1vmxgqy">No arbitrary parameters (unlike ORA, no fiddle parameters for &quot;foreground genes&quot;)</li> <li class="svelte-1vmxgqy">Bootstrapping approach always uses appropriate &quot;background set&quot; (unlike ORA)</li> <li class="svelte-1vmxgqy">Available as a R package and through this website (<i>GOAT online</i> in above menu)</li></ul> <span><i>*ORA: overrepresentation analysis (i.e. classical Fisher-exact or hypergeometric test)</i></span></div>`,ge,u,H,Ce="GOAT finds more significant genesets",he,L,qe="GOAT identifies more significant GO terms as compared to GSEA and ORA across 6 omics studies.",pe,N,w,we,ue,z,Ie=`Differential expression analysis results from each study (panels) were extracted from
			supplementary data and subjected to 5 approaches (colors) for identifying overrepresented GO
			terms. GO MF, GO CC and GO BP represent the respective Gene Ontology domains Molecular
			Functions, Cellular Components and Biological Processes. The x-axis shows the number of
			significant genesets identified by each method after Bonferroni adjustment at α = 0.05
			and adjustment for testing 3 sources of genesets.`,ve,v,F,De="GOAT p-values are accurate",me,S,Ee=`In simulations that generate 200 thousand random genelists and genesets, from small to large,
			we find that geneset p-values estimated by GOAT are accurate under the null hypothesis.`,fe,$,I,ke,xe,R,Ve=`The y-axis shows the observed geneset p-value for a given p-value threshold on the x-axis
			(both on -log10 scale, p-values are shown as-is without multiple testing correction). For
			example, selecting a threshold of p = 0.01 (2 on -log10 scale) on the x-axis shows the
			proportion of observed genesets with p ≤ 0.01 on the y-axis. Since randomized genesets were
			used in these simulations, expected values are on the diagonal (dashed line).`,ye,m,P,Me="Description of the GOAT algorithm",_e,K,He=`In brief, the Geneset Ordinal Association Test (GOAT) is a parameter-free permutation-based
			algorithm for geneset enrichment analysis. It is easy to use via the online webtool or R
			package, computationally efficient and the resulting geneset p-values are well calibrated
			under the null hypothesis and invariant to geneset size. Application to various real-world
			proteomics and gene expression studies demonstrates that GOAT consistently identifies more
			significant Gene Ontology terms as compared to alternative methods.`,be,J,U,ze,Ge,j,Fe=`<ol><li class="svelte-1vmxgqy">Required input are a list of genes and their respective test statistics (p-value /
					effectsize), and a list of genesets obtained from GO or alternative resources.</li> <li class="svelte-1vmxgqy">Test statistics from the genelist are transformed to gene scores by rank(-pvalue)^2 or
					rank(effectsize)^2 depending on user-input, i.e. smaller p-values translate to higher gene
					scores. The result is a skewed gene score distribution.</li> <li class="svelte-1vmxgqy">For each geneset size N (number of genes), bootstrapping procedures generate a null
					distribution of geneset scores. This yields a skewed-normal distribution for small
					genesets and converges to a normal distribution for large genesets.</li> <li class="svelte-1vmxgqy">Geneset significance is determined for each geneset by comparing its score (mean of
					respective gene scores) against a null distribution of the same size (N).</li></ol>`;return{c(){b=t("meta"),B=r(),i=t("nav"),p=t("div"),f=t("a"),y=t("img"),ae=r(),A=t("a"),ne=ee("Home"),ie=r(),O=t("a"),re=ee("GOAT online"),le=r(),T=t("a"),oe=ee("gene ID mapping"),de=r(),C=t("a"),ce=ee("Documentation"),Y=r(),q=t("div"),q.innerHTML=Oe,Z=r(),g=t("div"),M=t("div"),M.innerHTML=Te,ge=r(),u=t("div"),H=t("h2"),H.textContent=Ce,he=r(),L=t("div"),L.textContent=qe,pe=r(),N=t("div"),w=t("img"),ue=r(),z=t("p"),z.textContent=Ie,ve=r(),v=t("div"),F=t("h2"),F.textContent=De,me=r(),S=t("div"),S.textContent=Ee,fe=r(),$=t("div"),I=t("img"),xe=r(),R=t("p"),R.textContent=Ve,ye=r(),m=t("div"),P=t("h2"),P.textContent=Me,_e=r(),K=t("p"),K.textContent=He,be=r(),J=t("div"),U=t("img"),Ge=r(),j=t("div"),j.innerHTML=Fe,this.h()},l(c){const G=Je("svelte-1b313dc",document.head);b=s(G,"META",{name:!0,content:!0}),G.forEach(n),B=l(c),i=s(c,"NAV",{});var _=h(i);p=s(_,"DIV",{style:!0});var Re=h(p);f=s(Re,"A",{href:!0});var Pe=h(f);y=s(Pe,"IMG",{src:!0,width:!0,height:!0,alt:!0}),Pe.forEach(n),Re.forEach(n),ae=l(_),A=s(_,"A",{href:!0,style:!0});var je=h(A);ne=te(je,"Home"),je.forEach(n),ie=l(_),O=s(_,"A",{href:!0,style:!0});var Be=h(O);re=te(Be,"GOAT online"),Be.forEach(n),le=l(_),T=s(_,"A",{href:!0,style:!0});var Le=h(T);oe=te(Le,"gene ID mapping"),Le.forEach(n),de=l(_),C=s(_,"A",{href:!0,style:!0});var Ne=h(C);ce=te(Ne,"Documentation"),Ne.forEach(n),_.forEach(n),Y=l(c),q=s(c,"DIV",{class:!0,"data-svelte-h":!0}),x(q)!=="svelte-p77g5c"&&(q.innerHTML=Oe),Z=l(c),g=s(c,"DIV",{style:!0});var D=h(g);M=s(D,"DIV",{style:!0,"data-svelte-h":!0}),x(M)!=="svelte-1qljcv1"&&(M.innerHTML=Te),ge=l(D),u=s(D,"DIV",{style:!0});var E=h(u);H=s(E,"H2",{class:!0,"data-svelte-h":!0}),x(H)!=="svelte-10ytmrf"&&(H.textContent=Ce),he=l(E),L=s(E,"DIV",{"data-svelte-h":!0}),x(L)!=="svelte-t7brcz"&&(L.textContent=qe),pe=l(E),N=s(E,"DIV",{style:!0});var Se=h(N);w=s(Se,"IMG",{src:!0,width:!0,height:!0,alt:!0}),Se.forEach(n),ue=l(E),z=s(E,"P",{style:!0,"data-svelte-h":!0}),x(z)!=="svelte-8x3x4l"&&(z.textContent=Ie),E.forEach(n),ve=l(D),v=s(D,"DIV",{style:!0});var k=h(v);F=s(k,"H2",{class:!0,"data-svelte-h":!0}),x(F)!=="svelte-1cxefkq"&&(F.textContent=De),me=l(k),S=s(k,"DIV",{"data-svelte-h":!0}),x(S)!=="svelte-36zne0"&&(S.textContent=Ee),fe=l(k),$=s(k,"DIV",{style:!0});var $e=h($);I=s($e,"IMG",{src:!0,width:!0,height:!0,alt:!0}),$e.forEach(n),xe=l(k),R=s(k,"P",{style:!0,"data-svelte-h":!0}),x(R)!=="svelte-ogerj8"&&(R.textContent=Ve),k.forEach(n),ye=l(D),m=s(D,"DIV",{style:!0});var V=h(m);P=s(V,"H2",{class:!0,"data-svelte-h":!0}),x(P)!=="svelte-69l058"&&(P.textContent=Me),_e=l(V),K=s(V,"P",{"data-svelte-h":!0}),x(K)!=="svelte-150sgbv"&&(K.textContent=He),be=l(V),J=s(V,"DIV",{});var Ke=h(J);U=s(Ke,"IMG",{src:!0,alt:!0}),Ke.forEach(n),Ge=l(V),j=s(V,"DIV",{style:!0,"data-svelte-h":!0}),x(j)!=="svelte-rfve2k"&&(j.innerHTML=Fe),V.forEach(n),D.forEach(n),this.h()},h(){document.title="GOAT: Geneset Ordinal Association Test",a(b,"name","description"),a(b,"content","Geneset enrichment analysis for Geneset Ontology (GO) or KEGG pathways using the GOAT algorithm webtool. Online data analysis for your preranked genelist from e.g. proteomics or bulk/scRNAseq gene expression studies"),se(y.src,W=o[0]+"android-chrome-192x192.png")||a(y,"src",W),a(y,"width","40"),a(y,"height","40"),a(y,"alt","GOAT"),a(f,"href",o[0]),Q(f,"active",o[1]==="home"),d(p,"padding","4px"),d(p,"margin-left","20px"),a(A,"href",o[0]),d(A,"margin-left","5px"),Q(A,"active",o[1]==="home"),a(O,"href",o[0]+"goat"),d(O,"margin-left","40px"),Q(O,"active",o[1]==="goat"),a(T,"href",o[0]+"genemap"),d(T,"margin-left","40px"),Q(T,"active",o[1]==="genemap"),a(C,"href",o[0]+"docs"),d(C,"margin-left","40px"),Q(C,"active",o[1]==="docs"),a(q,"class","divTitle svelte-1vmxgqy"),d(M,"padding-top","50px"),a(H,"class","svelte-1vmxgqy"),se(w.src,we=o[0]+"barplot_signif_count.svg")||a(w,"src",we),a(w,"width","800"),a(w,"height","450"),a(w,"alt","barplot significant geneset counts"),d(N,"margin-top","25px"),d(z,"padding","0px 50px 0px 20px"),d(u,"padding-top","50px"),a(F,"class","svelte-1vmxgqy"),se(I.src,ke=o[0]+"null_simulations.svg")||a(I,"src",ke),a(I,"width","800"),a(I,"height","225"),a(I,"alt","GOAT p-values are accurate under the null hypothesis"),d($,"margin-top","25px"),d(R,"padding","0px 50px 0px 20px"),d(v,"padding-top","50px"),a(P,"class","svelte-1vmxgqy"),se(U.src,ze=o[0]+"goat_algorithm.png")||a(U,"src",ze),a(U,"alt","GOAT algorithm"),d(j,"padding","0px 50px 0px 10px"),d(m,"margin-top","100px"),d(g,"background-color","white"),d(g,"padding","0px 50px 0px 50px")},m(c,G){e(document.head,b),X(c,B,G),X(c,i,G),e(i,p),e(p,f),e(f,y),e(i,ae),e(i,A),e(A,ne),e(i,ie),e(i,O),e(O,re),e(i,le),e(i,T),e(T,oe),e(i,de),e(i,C),e(C,ce),X(c,Y,G),X(c,q,G),X(c,Z,G),X(c,g,G),e(g,M),e(g,ge),e(g,u),e(u,H),e(u,he),e(u,L),e(u,pe),e(u,N),e(N,w),e(u,ue),e(u,z),e(g,ve),e(g,v),e(v,F),e(v,me),e(v,S),e(v,fe),e(v,$),e($,I),e(v,xe),e(v,R),e(g,ye),e(g,m),e(m,P),e(m,_e),e(m,K),e(m,be),e(m,J),e(J,U),e(m,Ge),e(m,j)},p:Ae,i:Ae,o:Ae,d(c){c&&(n(B),n(i),n(Y),n(q),n(Z),n(g)),n(b)}}}function tt(o,b,B){let i;Qe(o,Ze,W=>B(2,i=W));const p=Ue&&Ue+"/"||!!i&&i.url.pathname==="/"&&"/"||"/goat/",f=i&&i.url?i.url.pathname.replace(".html",""):"/",y=f===p+"/goat"&&"goat"||f===p+"/docs"&&"docs"||f===p+"/genemap"&&"genemap"||"home";return[p,y]}class lt extends We{constructor(b){super(),Ye(this,b,tt,et,Xe,{})}}export{lt as component};