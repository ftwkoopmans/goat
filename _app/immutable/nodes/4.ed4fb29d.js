import{s as mn,p as Jt,E as bn,f as d,a as y,l as P,A as vn,g as c,d as n,c as w,h as A,B as V,m as W,j as _,D as Ie,k as B,F as Xe,x as s,i as b,G,H as Qt,I as gn,z as $e,n as Ce,y as fn,e as Xt}from"../chunks/scheduler.a30ae394.js";import{S as _n,i as yn,f as $t,b as en,d as tn,m as nn,a as ln,t as sn,e as an}from"../chunks/index.49b3fccf.js";import{p as wn,S as on,g as In,f as Cn,a as kn,d as Tn,c as xn,F as An,l as En}from"../chunks/stores.6ea4631d.js";/* empty css                    */import{p as Rn}from"../chunks/stores.b6b25a96.js";import{e as ce,w as rn}from"../chunks/singletons.3d14f100.js";function Gn(e){const t=e.length,r=e[0].map(a=>a.toLowerCase()),f=r.indexOf("symbol");let l;const p=f+1,i=[];if(f===-1)throw new Error("input genelist table must contain headers on the first line. Cannot find required header for column; 'symbol' (should contain character/string values)");if(r.filter(a=>a==="symbol").length>1)throw new Error("input genelist table is invalid because it contains multiple columns with the name 'symbol' (remove redundant columns)");for(let a=1;a<t;a++)if(e[a].length===0)i.push("");else{if(l=e[a],l.length<p)throw new Error("line "+(a+1)+" is invalid, it contains fewer values/columns than expected");i.push(l[f].toUpperCase())}return i}function Dn(e,t,r){const f=new RegExp("[ ,;]"),l=[...new Set(e.filter(a=>a.length>=2&&f.test(a)===!1))],p=[],i=e.length;for(let a=0;a<i;a++)if(e[a]==="")p.push([]);else{const h=[...new Set(e[a].split(/[ ,;]+/))].filter(m=>m.length>=2),T=h.length;let I;if(T===1)I=t[h[0]];else if(T>1&&r!=="skipAmbiguous")if(r==="retainFirstSucces")for(let m=0;m<T&&I===void 0;m++)I=t[h[m]];else if(r==="retainFirstSuccesNovel")for(let m=0;m<T&&I===void 0;m++)l.indexOf(h[m])===-1&&(I=t[h[m]]);else r==="retainFirst"&&(I=t[h[0]]);p.push(I===void 0?"":I)}return p}async function Sn(e){const t=await wn(e);if(Array.isArray(t)===!1||t.length<2)throw new Error("input genelist table must not be empty");return t}async function Hn(e,t,r,f){let l=[...e];const p=l.length,i=Dn(t,r.lookupGene,f);let a=[],h;for(let u=0;u<l[0].length;u++)h=l[0][u],h&&["gene","hgnc_id","hgnc_symbol"].indexOf(h.toLowerCase())!==-1&&a.push(u);const T=a.length;if(T>0){a=a.reverse();for(let u=0;u<p;u++)for(let R=0;R<T;R++)l[u].splice(a[R],1)}l[0]=l[0].concat(["gene","hgnc_id","hgnc_symbol"]);let I,m;for(let u=1;u<p;u++)I=i[u-1],I===""?l[u]=l[u].concat(["","",""]):(m=r.lookupMeta[I],m===void 0?l[u]=l[u].concat(["","",""]):l[u]=l[u].concat([m.entrez_id,I,m.hgnc_symbol]));return l}const{document:Et}=In;function Vn(e){let t,r="status:",f,l,p,i,a="no genelist loaded yet";return{c(){t=d("b"),t.textContent=r,f=y(),l=d("br"),p=y(),i=d("div"),i.textContent=a,this.h()},l(h){t=c(h,"B",{"data-svelte-h":!0}),V(t)!=="svelte-jrlymk"&&(t.textContent=r),f=w(h),l=c(h,"BR",{}),p=w(h),i=c(h,"DIV",{class:!0,"data-svelte-h":!0}),V(i)!=="svelte-1x7jbfu"&&(i.textContent=a),this.h()},h(){_(i,"class","alert")},m(h,T){b(h,t,T),b(h,f,T),b(h,l,T),b(h,p,T),b(h,i,T)},p:fn,d(h){h&&(n(t),n(f),n(l),n(p),n(i))}}}function Ln(e){let t,r="filename:",f,l=e[1].filename+"",p,i,a,h,T="number of rows with a gene symbol:",I,m=e[1].countGenes+"",u,R;return{c(){t=d("b"),t.textContent=r,f=y(),p=P(l),i=d("br"),a=y(),h=d("b"),h.textContent=T,I=y(),u=P(m),R=d("br")},l(v){t=c(v,"B",{"data-svelte-h":!0}),V(t)!=="svelte-nlrd5p"&&(t.textContent=r),f=w(v),p=W(v,l),i=c(v,"BR",{}),a=w(v),h=c(v,"B",{"data-svelte-h":!0}),V(h)!=="svelte-1vlyilb"&&(h.textContent=T),I=w(v),u=W(v,m),R=c(v,"BR",{})},m(v,E){b(v,t,E),b(v,f,E),b(v,p,E),b(v,i,E),b(v,a,E),b(v,h,E),b(v,I,E),b(v,u,E),b(v,R,E)},p(v,E){E[0]&2&&l!==(l=v[1].filename+"")&&Ce(p,l),E[0]&2&&m!==(m=v[1].countGenes+"")&&Ce(u,m)},d(v){v&&(n(t),n(f),n(p),n(i),n(a),n(h),n(I),n(u),n(R))}}}function Bn(e){let t,r="error while loading your genelist:",f,l,p;return{c(){t=d("b"),t.textContent=r,f=y(),l=d("div"),p=P(e[3]),this.h()},l(i){t=c(i,"B",{"data-svelte-h":!0}),V(t)!=="svelte-rtmaud"&&(t.textContent=r),f=w(i),l=c(i,"DIV",{class:!0});var a=A(l);p=W(a,e[3]),a.forEach(n),this.h()},h(){_(l,"class","alert")},m(i,a){b(i,t,a),b(i,f,a),b(i,l,a),s(l,p)},p(i,a){a[0]&8&&Ce(p,i[3])},d(i){i&&(n(t),n(f),n(l))}}}function un(e){let t,r="error while comparing your data to the lookup tables:",f,l,p;return{c(){t=d("b"),t.textContent=r,f=y(),l=d("div"),p=P(e[4]),this.h()},l(i){t=c(i,"B",{"data-svelte-h":!0}),V(t)!=="svelte-133bokb"&&(t.textContent=r),f=w(i),l=c(i,"DIV",{class:!0});var a=A(l);p=W(a,e[4]),a.forEach(n),this.h()},h(){_(l,"class","alert")},m(i,a){b(i,t,a),b(i,f,a),b(i,l,a),s(l,p)},p(i,a){a[0]&16&&Ce(p,i[4])},d(i){i&&(n(t),n(f),n(l))}}}function Nn(e){let t,r="No results yet; load your genelist and press START...";return{c(){t=d("div"),t.textContent=r},l(f){t=c(f,"DIV",{"data-svelte-h":!0}),V(t)!=="svelte-fomj7a"&&(t.textContent=r)},m(f,l){b(f,t,l)},p:fn,d(f){f&&n(t)}}}function Mn(e){let t,r=e[1].countSuccess+"",f,l,p=e[1].countGenes+"",i,a,h,T,I,m,u,R="download results",v,E,ee,j,S=!!e[5]&&dn(e);return{c(){t=d("div"),f=P(r),l=P(" / "),i=P(p),a=P(` genes in your table were mapped to Human gene\r
				identifiers.`),h=y(),T=d("br"),I=y(),m=d("div"),u=d("button"),u.textContent=R,v=y(),S&&S.c(),E=Xt(),this.h()},l(C){t=c(C,"DIV",{});var x=A(t);f=W(x,r),l=W(x," / "),i=W(x,p),a=W(x,` genes in your table were mapped to Human gene\r
				identifiers.`),x.forEach(n),h=w(C),T=c(C,"BR",{}),I=w(C),m=c(C,"DIV",{});var he=A(m);u=c(he,"BUTTON",{class:!0,"data-svelte-h":!0}),V(u)!=="svelte-c9mrqp"&&(u.textContent=R),he.forEach(n),v=w(C),S&&S.l(C),E=Xt(),this.h()},h(){_(u,"class","btn btn-action")},m(C,x){b(C,t,x),s(t,f),s(t,l),s(t,i),s(t,a),b(C,h,x),b(C,T,x),b(C,I,x),b(C,m,x),s(m,u),b(C,v,x),S&&S.m(C,x),b(C,E,x),ee||(j=G(u,"click",e[15]),ee=!0)},p(C,x){x[0]&2&&r!==(r=C[1].countSuccess+"")&&Ce(f,r),x[0]&2&&p!==(p=C[1].countGenes+"")&&Ce(i,p),C[5]?S?S.p(C,x):(S=dn(C),S.c(),S.m(E.parentNode,E)):S&&(S.d(1),S=null)},d(C){C&&(n(t),n(h),n(T),n(I),n(m),n(v),n(E)),S&&S.d(C),ee=!1,j()}}}function dn(e){let t,r="error while preparing data for download:",f,l,p;return{c(){t=d("b"),t.textContent=r,f=y(),l=d("div"),p=P(e[5]),this.h()},l(i){t=c(i,"B",{"data-svelte-h":!0}),V(t)!=="svelte-ncd7xk"&&(t.textContent=r),f=w(i),l=c(i,"DIV",{class:!0});var a=A(l);p=W(a,e[5]),a.forEach(n),this.h()},h(){_(l,"class","alert")},m(i,a){b(i,t,a),b(i,f,a),b(i,l,a),s(l,p)},p(i,a){a[0]&32&&Ce(p,i[5])},d(i){i&&(n(t),n(f),n(l))}}}function Fn(e){let t,r,f,l,p,i=`<img src="${ce+"/android-chrome-192x192.png"}" width="40" height="40" alt="GOAT"/>`,a,h,T="Home",I,m,u="GOAT online",R,v,E="gene ID mapping",ee,j,S="Documentation",C,x,he=`<h1>Gene ID mapping</h1> <div style="margin: 0px 25px 0px 25px; text-align: center;"><h3>With this tool you can map gene symbols in an Excel/CSV/TSV table to Human gene identifiers.
			<br/>
			The output will contain columns with gene IDs and the respective official gene symbol. It can be
			directly used as input for the GOAT online tool (it&#39;ll only use rows where mapping was successful).</h3> <img src="beta.png" width="200px" alt="BETA VERSION"/></div>`,ke,D,ae,qe="<h2>Your input data</h2>",Ne,me,be,oe,J,Q,Me,N,K,g,H="drag&drop a genelist file or click to open a file dialog",te,fe,et,Fe,tt,re,nt,lt,Te,Rt=`<div class="helpSettings">Importantly, your input genelist / gene table must be prepared in a format that is
						compatible with this tool.
						<ul><li class="svelte-1b7cd8l">File format: either CSV, TSV or Excel (.xlsx file, data on the first sheet)</li> <li class="svelte-1b7cd8l">Required column: <b>symbol</b> (column name must match exactly)</li></ul>
						The documentation below shows an example table that one might use as input.</div>`,st,ve,Oe,L,ge,Gt="<h2>Settings</h2>",it,_e,Dt=`Rows where the gene symbol column contains a delimiter (i.e. semicolon, comma or whitespace) are\r
		assumed to refer to multiple genes. If there are more than 2 unique gene symbols (on a row), how\r
		should these be mapped to human gene identifiers?`,at,xe,Ae,q,ot,rt,Ee,Re,O,ut,dt,Ge,De,z,ct,ft,Se,He,U,pt,ht,mt,bt,ne,le,vt,ze,gt,ue,_t,yt,Ue,se,ye,St="<h2>Results</h2>",wt,we,Pe,ie,Ht=`<div class="divH2"><h2>Settings documentation</h2></div> <p>To illustrate the problem of ambiguous genes/symbols and solutions offered by above options;</p> <table class="tbl-bordered"><thead><tr><th>symbol</th><th>effectsize</th><th class="note svelte-1b7cd8l">note</th></tr></thead> <tbody><tr><td>GRIA1</td><td>1.0</td> <td class="note svelte-1b7cd8l">protein group maps to exactly 1 gene</td></tr> <tr><td>GRIA2</td><td>1.0</td> <td class="note svelte-1b7cd8l">protein group maps to exactly 1 gene</td></tr> <tr><td>GRIA1;GRIA2</td><td>1.5</td> <td class="note svelte-1b7cd8l">ambiguous, one might want to use only the first entry (&#39;leading&#39; gene)</td></tr> <tr><td>GRIA3;GRIA4</td><td>1.5</td> <td class="note svelte-1b7cd8l">ambiguous, but this row contributes a new gene (GRIA3)</td></tr> <tr><td>tr|A8K0K0|A8K0K0_HUMAN;GRIA2;GRIA3</td><td>2.0</td> <td class="note svelte-1b7cd8l">ambiguous, but the first entry has no gene symbol only an accession</td></tr></tbody></table> <div><ol><li class="svelte-1b7cd8l">With option 1, <i>&quot;try to map only the first gene symbol&quot;</i><br/>
				the output will contain a gene ID for all rows except the last (respectively, GRIA1, GRIA2, GRIA1,
				GRIA3, -).</li> <li class="svelte-1b7cd8l">With option 2, <i>&quot;use the first gene symbol that can be successfully mapped&quot;</i><br/>
				all rows will be mapped to a gene ID (respectively, GRIA1, GRIA2, GRIA1, GRIA3, GRIA2).</li> <li class="svelte-1b7cd8l">With option 3, <i>&quot;analogous to above, but skip if there exist a non-ambiguous row that contains the first
					(successfully mapped) gene symbol&quot;</i><br/>
				rows 3 and 5 are skipped because there exist unambiguous entries for GRIA1 and GRIA2 (respectively,
				GRIA1, GRIA2, -, GRIA3, -). This approach favors rows that are unambiguous and supplements this
				only with ambiguous rows that contribute new information (genes).</li> <li class="svelte-1b7cd8l">With option 4, <i>&quot;skip rows with ambiguous symbols altogether&quot;</i><br/>
				only the first 2 rows are mapped (respectively, GRIA1, GRIA2, -, -, -).</li></ol></div> <p>When chosing the option most appropriate for your dataset, keep in mind that the geneset
		analysis in GOAT online will retain only 1 row per unique gene. If multiple rows/entries are
		available for a gene, the one with the lowest/best p-value is retained. If there are no gene
		p-values in your data/table, the best absolute effectsize is retained (across multiple entries
		for the same gene).</p>`,We,de,Vt=`<div class="divH2"><h2>What data is used / how is gene ID mapping done?</h2></div> <div><p>We created a lookup table using official gene symbols and aliases/synonyms based on
			information provided by <a href="https://www.genenames.org" target="_blank" rel="nofollow">HGNC</a>, it is stored on this webserver. Any synonyms that are listed by HGNC as entries for
			multiple genes are considered ambiguous and are discarded.</p> <p>This tool downloads the lookup table (mapping from symbols to gene IDs) to your computer and
			then proceeds with matching between your input table and HGNC gene information. So your
			table/data never leaves your computer and remains private at all times.</p></div>`,Ve,It,Ct,Lt;function pn(o){e[17](o)}let Bt={colour:"black",size:"2.5em"};e[6]!==void 0&&(Bt.show=e[6]),re=new on({props:Bt}),Jt.push(()=>$t(re,"show",pn));function Nt(o,k){return o[3]?Bn:typeof o[1]=="object"&&o[1].ok===!0?Ln:Vn}let Ke=Nt(e),X=Ke(e);function hn(o){e[28](o)}let Mt={colour:"black",size:"1.5em"};e[7]!==void 0&&(Mt.show=e[7]),ue=new on({props:Mt}),Jt.push(()=>$t(ue,"show",hn));let M=!!e[4]&&un(e);function Ft(o,k){return o[7]===!1&&typeof o[1]=="object"&&o[1].done===!0?Mn:Nn}let Ye=Ft(e),$=Ye(e);return It=bn(e[24][0]),{c(){t=d("meta"),r=y(),f=d("nav"),l=d("div"),p=d("a"),p.innerHTML=i,a=y(),h=d("a"),h.textContent=T,I=y(),m=d("a"),m.textContent=u,R=y(),v=d("a"),v.textContent=E,ee=y(),j=d("a"),j.textContent=S,C=y(),x=d("div"),x.innerHTML=he,ke=y(),D=d("div"),ae=d("div"),ae.innerHTML=qe,Ne=y(),me=d("table"),be=d("tbody"),oe=d("tr"),J=d("td"),Q=d("input"),Me=y(),N=d("div"),K=d("p"),g=d("i"),g.textContent=H,te=y(),fe=d("br"),et=y(),Fe=d("i"),tt=y(),en(re.$$.fragment),lt=y(),Te=d("td"),Te.innerHTML=Rt,st=y(),ve=d("div"),X.c(),Oe=y(),L=d("div"),ge=d("div"),ge.innerHTML=Gt,it=y(),_e=d("div"),_e.textContent=Dt,at=y(),xe=d("div"),Ae=d("label"),q=d("input"),ot=P(`\r
\r
			try to map only the first gene symbol`),rt=y(),Ee=d("div"),Re=d("label"),O=d("input"),ut=P(`\r
\r
			use the first gene symbol that can be successfully mapped`),dt=y(),Ge=d("div"),De=d("label"),z=d("input"),ct=P(`\r
\r
			analogous to above, but skip if there exist a non-ambiguous row that contains the first\r
			(successfully mapped) gene symbol`),ft=y(),Se=d("div"),He=d("label"),U=d("input"),pt=P(`\r
\r
			skip rows with ambiguous symbols altogether`),ht=y(),mt=d("br"),bt=y(),ne=d("div"),le=d("button"),vt=P("START"),gt=y(),en(ue.$$.fragment),yt=y(),M&&M.c(),Ue=y(),se=d("div"),ye=d("div"),ye.innerHTML=St,wt=y(),we=d("div"),$.c(),Pe=y(),ie=d("div"),ie.innerHTML=Ht,We=y(),de=d("div"),de.innerHTML=Vt,this.h()},l(o){const k=vn("svelte-1bjj788",Et.head);t=c(k,"META",{name:!0,content:!0}),k.forEach(n),r=w(o),f=c(o,"NAV",{});var Y=A(f);l=c(Y,"DIV",{style:!0});var je=A(l);p=c(je,"A",{href:!0,"data-svelte-h":!0}),V(p)!=="svelte-6nvu5r"&&(p.innerHTML=i),je.forEach(n),a=w(Y),h=c(Y,"A",{href:!0,style:!0,"data-svelte-h":!0}),V(h)!=="svelte-xsa54j"&&(h.textContent=T),I=w(Y),m=c(Y,"A",{href:!0,style:!0,"data-svelte-h":!0}),V(m)!=="svelte-o349n8"&&(m.textContent=u),R=w(Y),v=c(Y,"A",{href:!0,style:!0,"data-svelte-h":!0}),V(v)!=="svelte-sn5cge"&&(v.textContent=E),ee=w(Y),j=c(Y,"A",{href:!0,style:!0,"data-svelte-h":!0}),V(j)!=="svelte-1dikf5g"&&(j.textContent=S),Y.forEach(n),C=w(o),x=c(o,"DIV",{class:!0,"data-svelte-h":!0}),V(x)!=="svelte-15c5tkv"&&(x.innerHTML=he),ke=w(o),D=c(o,"DIV",{class:!0,style:!0,tabindex:!0,role:!0,"aria-pressed":!0});var Le=A(D);ae=c(Le,"DIV",{class:!0,"data-svelte-h":!0}),V(ae)!=="svelte-kki37h"&&(ae.innerHTML=qe),Ne=w(Le),me=c(Le,"TABLE",{});var jt=A(me);be=c(jt,"TBODY",{});var qt=A(be);oe=c(qt,"TR",{});var Ze=A(oe);J=c(Ze,"TD",{style:!0});var Je=A(J);Q=c(Je,"INPUT",{id:!0,type:!0,accept:!0,style:!0}),Me=w(Je),N=c(Je,"DIV",{id:!0,tabindex:!0,role:!0,"aria-pressed":!0});var Ot=A(N);K=c(Ot,"P",{});var pe=A(K);g=c(pe,"I",{"data-svelte-h":!0}),V(g)!=="svelte-uobxa9"&&(g.textContent=H),te=w(pe),fe=c(pe,"BR",{}),et=w(pe),Fe=c(pe,"I",{class:!0}),A(Fe).forEach(n),tt=w(pe),tn(re.$$.fragment,pe),pe.forEach(n),Ot.forEach(n),Je.forEach(n),lt=w(Ze),Te=c(Ze,"TD",{"data-svelte-h":!0}),V(Te)!=="svelte-i76nm3"&&(Te.innerHTML=Rt),Ze.forEach(n),qt.forEach(n),jt.forEach(n),st=w(Le),ve=c(Le,"DIV",{style:!0});var zt=A(ve);X.l(zt),zt.forEach(n),Le.forEach(n),Oe=w(o),L=c(o,"DIV",{class:!0,style:!0});var F=A(L);ge=c(F,"DIV",{class:!0,"data-svelte-h":!0}),V(ge)!=="svelte-1etyc0x"&&(ge.innerHTML=Gt),it=w(F),_e=c(F,"DIV",{style:!0,"data-svelte-h":!0}),V(_e)!=="svelte-atx6a9"&&(_e.textContent=Dt),at=w(F),xe=c(F,"DIV",{style:!0});var Ut=A(xe);Ae=c(Ut,"LABEL",{});var kt=A(Ae);q=c(kt,"INPUT",{type:!0,name:!0}),ot=W(kt,`\r
\r
			try to map only the first gene symbol`),kt.forEach(n),Ut.forEach(n),rt=w(F),Ee=c(F,"DIV",{style:!0});var Pt=A(Ee);Re=c(Pt,"LABEL",{});var Tt=A(Re);O=c(Tt,"INPUT",{type:!0,name:!0}),ut=W(Tt,`\r
\r
			use the first gene symbol that can be successfully mapped`),Tt.forEach(n),Pt.forEach(n),dt=w(F),Ge=c(F,"DIV",{style:!0});var Wt=A(Ge);De=c(Wt,"LABEL",{});var xt=A(De);z=c(xt,"INPUT",{type:!0,name:!0}),ct=W(xt,`\r
\r
			analogous to above, but skip if there exist a non-ambiguous row that contains the first\r
			(successfully mapped) gene symbol`),xt.forEach(n),Wt.forEach(n),ft=w(F),Se=c(F,"DIV",{style:!0});var Kt=A(Se);He=c(Kt,"LABEL",{});var At=A(He);U=c(At,"INPUT",{type:!0,name:!0}),pt=W(At,`\r
\r
			skip rows with ambiguous symbols altogether`),At.forEach(n),Kt.forEach(n),ht=w(F),mt=c(F,"BR",{}),bt=w(F),ne=c(F,"DIV",{});var Be=A(ne);le=c(Be,"BUTTON",{class:!0,style:!0});var Yt=A(le);vt=W(Yt,"START"),Yt.forEach(n),gt=w(Be),tn(ue.$$.fragment,Be),yt=w(Be),M&&M.l(Be),Be.forEach(n),F.forEach(n),Ue=w(o),se=c(o,"DIV",{class:!0,style:!0});var Qe=A(se);ye=c(Qe,"DIV",{class:!0,"data-svelte-h":!0}),V(ye)!=="svelte-1lm4it4"&&(ye.innerHTML=St),wt=w(Qe),we=c(Qe,"DIV",{style:!0});var Zt=A(we);$.l(Zt),Zt.forEach(n),Qe.forEach(n),Pe=w(o),ie=c(o,"DIV",{class:!0,style:!0,"data-svelte-h":!0}),V(ie)!=="svelte-1ydgkiu"&&(ie.innerHTML=Ht),We=w(o),de=c(o,"DIV",{class:!0,style:!0,"data-svelte-h":!0}),V(de)!=="svelte-1f7z8uw"&&(de.innerHTML=Vt),this.h()},h(){Et.title="Map gene symbols to human gene identifiers",_(t,"name","description"),_(t,"content","Map gene symbols in your input table to human gene identifiers (NCBI Entrez and HGNC)"),_(p,"href",ce+"/"),Ie(p,"active",e[8]==="home"),B(l,"padding","4px"),B(l,"margin-left","20px"),_(h,"href",ce+"/"),B(h,"margin-left","5px"),Ie(h,"active",e[8]==="home"),_(m,"href",ce+"/goat"),B(m,"margin-left","40px"),Ie(m,"active",e[8]==="goat"),_(v,"href",ce+"/genemap"),B(v,"margin-left","40px"),Ie(v,"active",e[8]==="genemap"),_(j,"href",ce+"/docs"),B(j,"margin-left","40px"),Ie(j,"active",e[8]==="docs"),_(x,"class","divH1"),_(ae,"class","divH2"),_(Q,"id","genelistFileDialog"),_(Q,"type","file"),_(Q,"accept",".csv,.tsv,.xlsx"),Q.multiple="false",B(Q,"display","none"),_(Fe,"class","dropicon"),_(N,"id","genelist_dropzone"),_(N,"tabindex","0"),_(N,"role","button"),_(N,"aria-pressed","false"),Ie(N,"is-dragover",e[0]===!0),B(J,"vertical-align","top"),B(J,"width","300px"),_(ve,"style",""),_(D,"class","pnlResults"),B(D,"margin-bottom","50px"),_(D,"tabindex","-1"),_(D,"role","button"),_(D,"aria-pressed","false"),_(ge,"class","divH2"),B(_e,"margin","10px 0px 10px 0px"),_(q,"type","radio"),_(q,"name","inputSymbolRetain"),q.__value="retainFirst",Xe(q,q.__value),B(xe,"margin","5px"),_(O,"type","radio"),_(O,"name","inputSymbolRetain"),O.__value="retainFirstSucces",Xe(O,O.__value),B(Ee,"margin","5px"),_(z,"type","radio"),_(z,"name","inputSymbolRetain"),z.__value="retainFirstSuccesNovel",Xe(z,z.__value),B(Ge,"margin","5px"),_(U,"type","radio"),_(U,"name","inputSymbolRetain"),U.__value="skipAmbiguous",Xe(U,U.__value),B(Se,"margin","5px"),_(le,"class","btn btn-action"),B(le,"margin-right","5px"),le.disabled=ze=e[7]===!0||e[6]===!0||!(typeof e[1]=="object"&&e[1].ok===!0),_(L,"class","pnlResults"),B(L,"margin-bottom","50px"),_(ye,"class","divH2"),B(we,"margin","10px 0px 10px 0px"),_(se,"class","pnlResults"),B(se,"margin-bottom","100px"),_(ie,"class","pnlResults"),B(ie,"margin-bottom","100px"),B(ie,"padding-right","50px"),_(de,"class","pnlResults"),B(de,"margin-bottom","100px"),It.p(q,O,z,U)},m(o,k){s(Et.head,t),b(o,r,k),b(o,f,k),s(f,l),s(l,p),s(f,a),s(f,h),s(f,I),s(f,m),s(f,R),s(f,v),s(f,ee),s(f,j),b(o,C,k),b(o,x,k),b(o,ke,k),b(o,D,k),s(D,ae),s(D,Ne),s(D,me),s(me,be),s(be,oe),s(oe,J),s(J,Q),s(J,Me),s(J,N),s(N,K),s(K,g),s(K,te),s(K,fe),s(K,et),s(K,Fe),s(K,tt),nn(re,K,null),s(oe,lt),s(oe,Te),s(D,st),s(D,ve),X.m(ve,null),b(o,Oe,k),b(o,L,k),s(L,ge),s(L,it),s(L,_e),s(L,at),s(L,xe),s(xe,Ae),s(Ae,q),q.checked=q.__value===e[2],s(Ae,ot),s(L,rt),s(L,Ee),s(Ee,Re),s(Re,O),O.checked=O.__value===e[2],s(Re,ut),s(L,dt),s(L,Ge),s(Ge,De),s(De,z),z.checked=z.__value===e[2],s(De,ct),s(L,ft),s(L,Se),s(Se,He),s(He,U),U.checked=U.__value===e[2],s(He,pt),s(L,ht),s(L,mt),s(L,bt),s(L,ne),s(ne,le),s(le,vt),s(ne,gt),nn(ue,ne,null),s(ne,yt),M&&M.m(ne,null),b(o,Ue,k),b(o,se,k),s(se,ye),s(se,wt),s(se,we),$.m(we,null),b(o,Pe,k),b(o,ie,k),b(o,We,k),b(o,de,k),Ve=!0,Ct||(Lt=[G(Q,"change",e[16]),G(N,"keydown",cn),G(N,"click",cn),G(N,"dragstart",Z),G(N,"dragover",e[18]),G(N,"dragenter",e[19]),G(N,"dragleave",e[20]),G(N,"dragend",e[21]),G(N,"drop",e[22]),G(D,"keydown",void 0),G(D,"dragstart",Z),G(D,"dragover",Z),G(D,"dragenter",Z),G(D,"dragleave",Z),G(D,"dragend",Z),G(D,"drop",Z),G(q,"change",e[13]),G(q,"change",e[23]),G(O,"change",e[13]),G(O,"change",e[25]),G(z,"change",e[13]),G(z,"change",e[26]),G(U,"change",e[13]),G(U,"change",e[27]),G(le,"click",e[14])],Ct=!0)},p(o,k){const Y={};!nt&&k[0]&64&&(nt=!0,Y.show=o[6],Qt(()=>nt=!1)),re.$set(Y),(!Ve||k[0]&1)&&Ie(N,"is-dragover",o[0]===!0),Ke===(Ke=Nt(o))&&X?X.p(o,k):(X.d(1),X=Ke(o),X&&(X.c(),X.m(ve,null))),k[0]&4&&(q.checked=q.__value===o[2]),k[0]&4&&(O.checked=O.__value===o[2]),k[0]&4&&(z.checked=z.__value===o[2]),k[0]&4&&(U.checked=U.__value===o[2]),(!Ve||k[0]&194&&ze!==(ze=o[7]===!0||o[6]===!0||!(typeof o[1]=="object"&&o[1].ok===!0)))&&(le.disabled=ze);const je={};!_t&&k[0]&128&&(_t=!0,je.show=o[7],Qt(()=>_t=!1)),ue.$set(je),o[4]?M?M.p(o,k):(M=un(o),M.c(),M.m(ne,null)):M&&(M.d(1),M=null),Ye===(Ye=Ft(o))&&$?$.p(o,k):($.d(1),$=Ye(o),$&&($.c(),$.m(we,null)))},i(o){Ve||(ln(re.$$.fragment,o),ln(ue.$$.fragment,o),Ve=!0)},o(o){sn(re.$$.fragment,o),sn(ue.$$.fragment,o),Ve=!1},d(o){o&&(n(r),n(f),n(C),n(x),n(ke),n(D),n(Oe),n(L),n(Ue),n(se),n(Pe),n(ie),n(We),n(de)),n(t),an(re),X.d(),an(ue),M&&M.d(),$.d(),It.r(),Ct=!1,gn(Lt)}}}function Z(e){e.preventDefault(),e.stopPropagation()}async function cn(e){document.getElementById("genelistFileDialog").click()}function jn(e,t,r){let f,l,p,i;$e(e,En,g=>r(29,f=g)),$e(e,Rn,g=>r(30,l=g));const a=l&&l.url?l.url.pathname.replace(".html",""):"/",h=a===ce+"/goat"&&"goat"||a===ce+"/docs"&&"docs"||a===ce+"/genemap"&&"genemap"||"home",T=rn(!1);$e(e,T,g=>r(6,p=g));const I=rn(!1);$e(e,I,g=>r(7,i=g));let m=!1,u,R="retainFirstSucces",v,E,ee;async function j(g){g.dataTransfer&&g.dataTransfer.files&&S(g.dataTransfer.files[0])}async function S(g){try{r(3,v=void 0),T.update(fe=>!0);const H=await Sn(g),te=Gn(H);r(1,u={genelist:H,symbols:te,filename:g.name,ok:!0,result:[],countGenes:H.length-1,countSuccess:0,done:!1,analytics:!1})}catch(H){console.error(H),r(3,v=H.message)}T.update(H=>!1)}function C(g){typeof u=="object"&&(r(1,u.result=[],u),r(1,u.countSuccess=0,u),r(1,u.done=!1,u))}async function x(g){if(u===void 0||typeof u!="object"||typeof u.done!="boolean"){r(4,E="nothing to do, genelist has not been loaded");return}try{r(4,E=void 0),I.update(H=>!0),r(1,u.done=!1,u),await Cn(),r(1,u.result=await Hn(u.genelist,u.symbols,f,R),u),r(1,u.countSuccess=u.result.reduce((H,te)=>H+(te[te.length-1]===""?0:1),-1),u),r(1,u.done=!0,u),u.analytics!==!0&&(kn("genemap"),r(1,u.analytics=!0,u))}catch(H){console.error(H),r(4,E=H.message)}I.update(H=>!1)}async function he(g){try{const H=Tn(),te=await xn(u.result);let fe=u.filename.replace(/\.[a-zA-Z]{3,4}$/,"")+"__GOAT-online_genemap_"+H+".xlsx";fe.length>180&&(fe="GOAT-online_genemap_"+H+".xlsx"),An.saveAs(te,fe)}catch(H){console.error(H),r(5,ee=H.message)}}const ke=[[]],D=g=>{g.target.files&&S(g.target.files[0])};function ae(g){p=g,T.set(p)}const qe=g=>{Z(g),r(0,m=!0)},Ne=g=>{Z(g),r(0,m=!0)},me=g=>{Z(g),r(0,m=!1)},be=g=>{Z(g),r(0,m=!1)},oe=g=>{Z(g),r(0,m=!1),j(g)};function J(){R=this.__value,r(2,R)}function Q(){R=this.__value,r(2,R)}function Me(){R=this.__value,r(2,R)}function N(){R=this.__value,r(2,R)}function K(g){i=g,I.set(i)}return[m,u,R,v,E,ee,p,i,h,T,I,j,S,C,x,he,D,ae,qe,Ne,me,be,oe,J,ke,Q,Me,N,K]}class Kn extends _n{constructor(t){super(),yn(this,t,jn,Fn,mn,{},null,[-1,-1])}}export{Kn as component};