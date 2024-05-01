import{s as qe,f as s,a as o,l as Y,A as Le,g as l,d as a,c,h as r,m as Z,B as j,j as t,C as ee,k as i,x as e,i as z,y as ye,z as Ce}from"../chunks/scheduler.cd38125d.js";import{S as Be,i as Se}from"../chunks/index.a6a2bf8f.js";/* empty css                    */import{p as ze}from"../chunks/stores.bacf73c7.js";function Fe(g){let f,B,n,G,y,A,_e,te,b,se,le,M,ae,ie,V,ne,re,H,oe,P,T,Ge=`<h1 class="svelte-12gnch1">GOAT online gene set enrichment analysis</h1> <div><h3 class="svelte-12gnch1">Identify enriched Gene Ontology terms within a list of gene effect sizes from e.g. proteomics
			or gene expression studies with the Gene set Ordinal Association Test (GOAT)</h3></div>`,K,h,S,Ae=`<h1 class="svelte-12gnch1">GOAT algorithm features:</h1> <div class="listChecks svelte-12gnch1"><ul class="svelte-12gnch1"><li class="svelte-12gnch1">Sensitive: more significant gene sets as compared to other methods</li> <li class="svelte-12gnch1">Accurate: gene set p-values are accurate under the null hypothesis</li> <li class="svelte-12gnch1">Fast: completes in seconds</li> <li class="svelte-12gnch1">No arbitrary parameters:
					<br/>unlike ORA, no fiddle parameters for &quot;&gt;k foreground genes&quot; or &quot;gene significance
					cutoff&quot;</li> <li class="svelte-12gnch1">Bootstrapping approach always uses the appropriate &quot;background set&quot;:
					<br/>avoids all pitfalls commonly observed for ORA, as outlined by
					<a href="https://doi.org/10.1371/journal.pcbi.1009935" target="_blank" rel="nofollow" tabindex="-1">Wijesooriya et al. (2022)</a></li> <li class="svelte-12gnch1">Available as an R package and online data analysis tool (<i>GOAT online</i> in the above menu)</li></ul> <span><i>*ORA: classical overrepresentation analysis (i.e. Fisher-exact or hypergeometric test)</i></span></div> <div><p>The GOAT algorithm is described in the preprint manuscript;
				<br/> <i>Koopmans, F. (2023). GOAT: efficient and robust identification of gene set enrichment.</i> <br/> <a href="https://doi.org/10.1101/2023.12.10.570979" target="_blank" rel="nofollow">https://doi.org/10.1101/2023.12.10.570979</a></p></div>`,ce,v,k,F,de,O,be,he,q,Te='<h2 class="svelte-12gnch1">Interactive data analysis tools</h2> <ul class="svelte-12gnch1"><li class="svelte-12gnch1">Easy to use</li> <li class="svelte-12gnch1">Gene set databases from GO, SynGO, or upload your own (e.g. KEGG)</li> <li class="svelte-12gnch1">Download results in convenient Excel format</li> <li class="svelte-12gnch1">Includes Methods text tailored to your analyses</li> <li class="svelte-12gnch1">Create customized publication-ready figures</li></ul>',ge,u,L,Oe='<h2 class="svelte-12gnch1">Customize your figures</h2> <ul class="svelte-12gnch1"><li class="svelte-12gnch1">Resize the figure, change text size, colors, etc.</li> <li class="svelte-12gnch1">Figures are live-updated as you change settings</li> <li class="svelte-12gnch1">High-resolution figures in SVG format</li></ul>',ve,x,N,ue,E,xe,pe,p,I,R,me,w,Ee,fe,C,Ie='<h2 class="svelte-12gnch1">A wide variety of data visualizations</h2> <ul class="svelte-12gnch1"><li class="svelte-12gnch1">Interactive data table</li> <li class="svelte-12gnch1">Barplots / Lollipop charts</li> <li class="svelte-12gnch1">GO relations as compact tree</li> <li class="svelte-12gnch1">Treemap of GO relations</li> <li class="svelte-12gnch1">Heatmap of gene set similarities</li></ul>';return{c(){f=s("meta"),B=o(),n=s("nav"),G=s("div"),y=s("a"),A=s("img"),te=o(),b=s("a"),se=Y("Home"),le=o(),M=s("a"),ae=Y("GOAT online"),ie=o(),V=s("a"),ne=Y("gene ID mapping"),re=o(),H=s("a"),oe=Y("Documentation"),P=o(),T=s("div"),T.innerHTML=Ge,K=o(),h=s("div"),S=s("div"),S.innerHTML=Ae,ce=o(),v=s("div"),k=s("div"),F=s("span"),de=o(),O=s("img"),he=o(),q=s("div"),q.innerHTML=Te,ge=o(),u=s("div"),L=s("div"),L.innerHTML=Oe,ve=o(),x=s("div"),N=s("span"),ue=o(),E=s("img"),pe=o(),p=s("div"),I=s("div"),R=s("span"),me=o(),w=s("img"),fe=o(),C=s("div"),C.innerHTML=Ie,this.h()},l(d){const _=Le("svelte-3dwj32",document.head);f=l(_,"META",{name:!0,content:!0}),_.forEach(a),B=c(d),n=l(d,"NAV",{});var m=r(n);G=l(m,"DIV",{style:!0});var we=r(G);y=l(we,"A",{href:!0,class:!0});var De=r(y);A=l(De,"IMG",{src:!0,width:!0,height:!0,alt:!0}),De.forEach(a),we.forEach(a),te=c(m),b=l(m,"A",{href:!0,style:!0,class:!0});var Me=r(b);se=Z(Me,"Home"),Me.forEach(a),le=c(m),M=l(m,"A",{href:!0,style:!0});var Ve=r(M);ae=Z(Ve,"GOAT online"),Ve.forEach(a),ie=c(m),V=l(m,"A",{href:!0,style:!0});var He=r(V);ne=Z(He,"gene ID mapping"),He.forEach(a),re=c(m),H=l(m,"A",{href:!0,style:!0});var ke=r(H);oe=Z(ke,"Documentation"),ke.forEach(a),m.forEach(a),P=c(d),T=l(d,"DIV",{class:!0,"data-svelte-h":!0}),j(T)!=="svelte-u8tpc"&&(T.innerHTML=Ge),K=c(d),h=l(d,"DIV",{class:!0});var D=r(h);S=l(D,"DIV",{"data-svelte-h":!0}),j(S)!=="svelte-1gfvj6w"&&(S.innerHTML=Ae),ce=c(D),v=l(D,"DIV",{style:!0});var $=r(v);k=l($,"DIV",{});var U=r(k);F=l(U,"SPAN",{class:!0}),r(F).forEach(a),de=c(U),O=l(U,"IMG",{src:!0,width:!0,style:!0,alt:!0}),U.forEach(a),he=c($),q=l($,"DIV",{style:!0,"data-svelte-h":!0}),j(q)!=="svelte-gusg0n"&&(q.innerHTML=Te),$.forEach(a),ge=c(D),u=l(D,"DIV",{style:!0});var W=r(u);L=l(W,"DIV",{style:!0,"data-svelte-h":!0}),j(L)!=="svelte-rd2rqq"&&(L.innerHTML=Oe),ve=c(W),x=l(W,"DIV",{style:!0});var X=r(x);N=l(X,"SPAN",{class:!0}),r(N).forEach(a),ue=c(X),E=l(X,"IMG",{src:!0,width:!0,style:!0,alt:!0}),X.forEach(a),W.forEach(a),pe=c(D),p=l(D,"DIV",{style:!0});var J=r(p);I=l(J,"DIV",{style:!0});var Q=r(I);R=l(Q,"SPAN",{class:!0}),r(R).forEach(a),me=c(Q),w=l(Q,"IMG",{src:!0,width:!0,style:!0,alt:!0}),Q.forEach(a),fe=c(J),C=l(J,"DIV",{style:!0,"data-svelte-h":!0}),j(C)!=="svelte-1r73910"&&(C.innerHTML=Ie),J.forEach(a),D.forEach(a),this.h()},h(){document.title="GOAT: Gene set Ordinal Association Test",t(f,"name","description"),t(f,"content","Gene set enrichment analysis for Gene Ontology (GO) or KEGG pathways using the GOAT algorithm web tool. Online data analysis for your preranked gene list from e.g. proteomics or bulk/scRNAseq gene expression studies"),ee(A.src,_e=g[0]+"android-chrome-192x192.png")||t(A,"src",_e),t(A,"width","40"),t(A,"height","40"),t(A,"alt","GOAT"),t(y,"href",g[0]),t(y,"class","active"),i(G,"padding","4px"),i(G,"margin-left","20px"),t(b,"href",g[0]),i(b,"margin-left","5px"),t(b,"class","active"),t(M,"href",g[0]+"goat"),i(M,"margin-left","40px"),t(V,"href",g[0]+"genemap"),i(V,"margin-left","40px"),t(H,"href",g[0]+"docs"),i(H,"margin-left","40px"),t(T,"class","divTitle svelte-12gnch1"),t(F,"class","valignBlock svelte-12gnch1"),ee(O.src,be=g[0]+"demo_treemap.svg")||t(O,"src",be),t(O,"width","400"),i(O,"vertical-align","middle"),t(O,"alt","Example treemap figure generated by GOAT online"),i(q,"padding","10px"),i(v,"margin-top","100px"),i(v,"display","grid"),i(v,"grid-template-columns","1fr 1fr"),i(L,"padding","10px"),t(N,"class","valignBlock svelte-12gnch1"),ee(E.src,xe=g[0]+"demo_barplot.png")||t(E,"src",xe),t(E,"width","450"),i(E,"vertical-align","middle"),t(E,"alt","Example treemap figure generated by GOAT online"),i(x,"text-align","center"),i(u,"margin-top","50px"),i(u,"display","grid"),i(u,"grid-template-columns","auto 1fr"),t(R,"class","valignBlock svelte-12gnch1"),ee(w.src,Ee=g[0]+"demo_heatmap.png")||t(w,"src",Ee),t(w,"width","300"),i(w,"vertical-align","middle"),t(w,"alt","Example heatmap figure generated by GOAT online"),t(I,"style",""),i(C,"padding","10px"),i(p,"margin-top","50px"),i(p,"display","grid"),i(p,"grid-template-columns","1fr auto"),t(h,"class","divMain svelte-12gnch1")},m(d,_){e(document.head,f),z(d,B,_),z(d,n,_),e(n,G),e(G,y),e(y,A),e(n,te),e(n,b),e(b,se),e(n,le),e(n,M),e(M,ae),e(n,ie),e(n,V),e(V,ne),e(n,re),e(n,H),e(H,oe),z(d,P,_),z(d,T,_),z(d,K,_),z(d,h,_),e(h,S),e(h,ce),e(h,v),e(v,k),e(k,F),e(k,de),e(k,O),e(v,he),e(v,q),e(h,ge),e(h,u),e(u,L),e(u,ve),e(u,x),e(x,N),e(x,ue),e(x,E),e(h,pe),e(h,p),e(p,I),e(I,R),e(I,me),e(I,w),e(p,fe),e(p,C)},p:ye,i:ye,o:ye,d(d){d&&(a(B),a(n),a(P),a(T),a(K),a(h)),a(f)}}}function Ne(g,f,B){let n;return Ce(g,ze,y=>B(1,n=y)),[!!n&&n.url.hostname.endsWith("github.io")&&"/goat/"||"/"]}class $e extends Be{constructor(f){super(),Se(this,f,Ne,Fe,qe,{})}}export{$e as component};
