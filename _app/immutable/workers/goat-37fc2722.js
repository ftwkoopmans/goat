(function(){"use strict";onmessage=e=>{const n=JSON.parse(e.data);try{k(n.gl,n.gs,n.nulldistr,n.scoreType,n.selectedMinSize,n.selectedMaxSize,n.selectedPvalueAdjustment,n.selectedPvalueThreshold)}catch(_){postMessage(_)}};async function k(e,n,_,f,a,r,t,s){const i=e.length,o=T(e,f);let c;const l=_.length;for(let g=0;g<l;g++)if(_[g].N===i){c=_[g];break}if(c===void 0)throw new Error("null distribution parameters for your genelist of length "+i+" are not available");const u=Q(n,e,+a,Math.min(+r,Math.ceil(i/2)),0);q(u,e,o,c);const p=B(u,t,+s);postMessage({scoretype_props:o,gs_filtered:u,stats:p,gl:e})}function T(e,n){function _(s,i){const o=s.length,c=o*o/1e3;let l=0;for(let u=0;u<o;u++)l=(o-u)*(o-u)/c,s[u][i]=l,l>s[u].sort_score&&(s[u].sort_score=l)}if(Array.isArray(e)===!1)throw new Error("genelist is expected to be an array (of objects, with properties like gene/symbol/pvalue/effectsize)");const f=e[0].hasOwnProperty("pvalue")&&typeof e[0].pvalue=="number",a=e[0].hasOwnProperty("effectsize")&&typeof e[0].effectsize=="number",r=e.length;for(let s=0;s<r;s++)e[s].sort_score=0,delete e[s].score_pvalue,delete e[s].score_effectsize_up,delete e[s].score_effectsize_down,delete e[s].score_effectsize_abs;let t=[];if(n==="effectsize"?(t=["effectsize_up","effectsize_down"],e.sort((s,i)=>i.gene-s.gene),f&&e.sort((s,i)=>s.pvalue-i.pvalue),e.sort((s,i)=>i.effectsize-s.effectsize),_(e,"score_effectsize_up"),e.sort((s,i)=>i.gene-s.gene),f&&e.sort((s,i)=>s.pvalue-i.pvalue),e.sort((s,i)=>s.effectsize-i.effectsize),_(e,"score_effectsize_down")):n==="effectsize_up"?(t=["effectsize_up"],e.sort((s,i)=>i.gene-s.gene),f&&e.sort((s,i)=>s.pvalue-i.pvalue),e.sort((s,i)=>i.effectsize-s.effectsize)):n==="effectsize_down"?(t=["effectsize_down"],e.sort((s,i)=>i.gene-s.gene),f&&e.sort((s,i)=>s.pvalue-i.pvalue),e.sort((s,i)=>s.effectsize-i.effectsize)):n==="effectsize_abs"?(t=["effectsize_abs"],e.sort((s,i)=>i.gene-s.gene),f&&e.sort((s,i)=>s.pvalue-i.pvalue),e.sort((s,i)=>Math.abs(i.effectsize)-Math.abs(s.effectsize)),_(e,"score_effectsize_abs")):n==="pvalue"&&(t=["pvalue"],e.sort((s,i)=>i.gene-s.gene),a&&e.sort((s,i)=>Math.abs(i.effectsize)-Math.abs(s.effectsize)),e.sort((s,i)=>s.pvalue-i.pvalue),_(e,"score_pvalue")),t.length===0)throw new Error("unknown scoretype parameter; "+n);return t}function D(e,n){const _=2/Math.sqrt(2*Math.PI),f=_*(n-1/n),a=Math.sqrt((1-_*_)*(n*n+1/(n*n))+2*(_*_)-1),r=e*a+f,t=2/(n+1/n);let s=n;r<0&&(s=1/s);let i=1;return r>=0?i=t*s*z(r/s,0,1,!1,!1):i=1-t*s*z(r/s,0,1,!0,!1),i>1&&(i=1),i}function I(e,n,_,f){return f<1.0001?z((e-n)/_,0,1,!1,!1):D((e-n)/_,f)}function P(e){const n=e.length,_=Array.from({length:n},(s,i)=>n-i),f=Array.from(e.keys()).sort((s,i)=>e[i]-e[s]),a=Array.from(f.keys()).sort((s,i)=>f[s]-f[i]);let r=new Float64Array(n),t=1;for(let s=0;s<n;s++)r[s]=n/_[s]*e[f[s]],t=t<r[s]?t:r[s],r[s]=t;for(let s=0;s<n;s++)e[s]=r[a[s]]>1?1:r[a[s]];return e}function F(e,n){const _=e.reduce((t,s)=>t+(typeof s=="number"&&!isNaN(s)),0),f=e.length;let a=0,r=new Float64Array(f);for(let t=0;t<f;t++)r[t]=NaN;if(_===0)return r;if(n==="bonferroni"){for(let t=0;t<f;t++)typeof e[t]=="number"&&!isNaN(e[t])&&(a=e[t]*_,r[t]=a>1?1:a);return r}if(n==="fdr"||n==="BH"){if(f===_)r=P(e);else{const t=e.filter(o=>typeof o=="number"&&!isNaN(o)),s=P(t);let i=0;for(let o=0;o<f;o++)typeof e[o]=="number"&&!isNaN(e[o])&&(r[o]=s[i],i+=1)}return r}console.error("unknown padjust method; "+n)}function q(e,n,_,f){const a=n.length,r={};for(let h=0;h<a;h++)r[n[h].gene]=n[h];const t=_.map(h=>"score_"+h),s=_.length,i=e.length;let o,c,l,u,p,g,m,N=new Array(s);for(let h=0;h<i;h++){let v=e[h].genesets.length;for(let d=0;d<v;d++){u=e[h].genesets[d].genes,p=u.length,e[h].genesets[d].score=NaN,e[h].genesets[d].score_type="",e[h].genesets[d].pvalue=NaN,e[h].genesets[d].pvalue_adjust=NaN,e[h].genesets[d].signif=void 0,N.fill(0),c=0,l="";for(let M=0;M<p;M++){if(m=r[u[M]],m===void 0)throw new Error("cannot find this geneset gene in genelist; "+u[M]);for(let y=0;y<s;y++)N[y]+=m[t[y]]}for(let M=0;M<s;M++)N[M]>c&&(c=N[M],l=_[M]);e[h].genesets[d].score=c/p,e[h].genesets[d].score_type=l,o=L(p,f),g=I(e[h].genesets[d].score,f.mu,o.sd,o.xi),g=Math.min(1,g*s),e[h].genesets[d].pvalue=g}}}function B(e,n,_){const f=e.length;let a=0;const r={signif_count:0,signif_max_ngenes:0,signif_max_ngenes_input:0,signif_max_ngenes_signif:0,signif_min_pvalue_adjust:1};for(let t=0;t<f;t++){e[t].stats={signif_count:0,signif_max_ngenes:0,signif_max_ngenes_input:0,signif_max_ngenes_signif:0,signif_min_pvalue_adjust:1};let s=e[t].genesets.length;if(s>0){if(e[t].genesets[0].hasOwnProperty("pvalue")===!1)throw new Error("genesets do not have a pvalue property, cannot compute adjusted p-values");let i=new Float64Array(s);for(let o=0;o<s;o++)i[o]=e[t].genesets[o].pvalue;i=F(i,n);for(let o=0;o<s;o++)a=i[o]*f,a=a>1?1:a,e[t].genesets[o].pvalue_adjust=a,e[t].genesets[o].signif=isNaN(a)?void 0:a<=_,e[t].genesets[o].signif===!0&&(e[t].stats.signif_count++,e[t].genesets[o].ngenes>e[t].stats.signif_max_ngenes&&(e[t].stats.signif_max_ngenes=e[t].genesets[o].ngenes),e[t].genesets[o].ngenes_input>e[t].stats.signif_max_ngenes_input&&(e[t].stats.signif_max_ngenes_input=e[t].genesets[o].ngenes_input),e[t].genesets[o].ngenes_signif>e[t].stats.signif_max_ngenes_signif&&(e[t].stats.signif_max_ngenes_signif=e[t].genesets[o].ngenes_signif),e[t].genesets[o].pvalue_adjust<e[t].stats.signif_min_pvalue_adjust&&(e[t].stats.signif_min_pvalue_adjust=e[t].genesets[o].pvalue_adjust));r.signif_count+=e[t].stats.signif_count,e[t].stats.signif_max_ngenes>r.signif_max_ngenes&&(r.signif_max_ngenes=e[t].stats.signif_max_ngenes),e[t].stats.signif_max_ngenes_input>r.signif_max_ngenes_input&&(r.signif_max_ngenes_input=e[t].stats.signif_max_ngenes_input),e[t].stats.signif_max_ngenes_signif>r.signif_max_ngenes_signif&&(r.signif_max_ngenes_signif=e[t].stats.signif_max_ngenes_signif),e[t].stats.signif_min_pvalue_adjust<r.signif_min_pvalue_adjust&&(r.signif_min_pvalue_adjust=e[t].stats.signif_min_pvalue_adjust)}}return r}function Q(e,n,_,f,a){if(Number.isInteger(_)===!1||_<10)throw new Error("minSize must be an integer number larger than 10");if(Number.isInteger(f)===!1||f<=_||f>Math.ceil(n.length/2))throw new Error("maxSize must be an integer number larger than minSize and smaller than half the genelist length");if((a===void 0||a===0)&&(a=9999999999999),Number.isInteger(a)===!1||a<f)throw new Error("maxSizeInput must be an integer number larger-equals than maxSize");const r=new Set(n.map(g=>g.gene)),t=new Set(n.filter(g=>g.signif===!0).map(g=>g.gene)),s=[],i=e.length;let o,c,l,u,p;for(let g=0;g<i;g++)if(o=e[g].genesets.length,o>0){p=[];for(let m=0;m<o;m++){if(c=e[g].genesets[m].genes,!Array.isArray(c))throw console.error(e[g].genesets[m]),new Error("invalid genesets object; no genes available for geneset at index "+m+" at source index "+g);c.length<a&&(l=c.filter(N=>r.has(N)),l.length>=_&&l.length<=f&&(u=l.filter(N=>t.has(N)),p.push({id:e[g].genesets[m].id,name:e[g].genesets[m].name,parents:e[g].genesets[m].parents,ngenes_input:c.length,ngenes:l.length,ngenes_signif:u.length,genes:l,genes_signif:u})))}s.push({source:e[g].source,source_version:e[g].source_version,genesets:p})}return s}function S(e,n){return n.sigma_0+n.sigma_1*e+n.sigma_2*Math.pow(e,2)+n.sigma_3*Math.pow(e,3)+n.sigma_4*Math.pow(e,4)+n.sigma_5*Math.pow(e,5)+n.sigma_6*Math.pow(e,6)+n.sigma_7*Math.pow(e,7)+n.sigma_8*Math.pow(e,8)+n.sigma_9*Math.pow(e,9)}function G(e,n){return n.xi_0+n.xi_1*e+n.xi_2*Math.pow(e,2)+n.xi_3*Math.pow(e,3)+n.xi_4*Math.pow(e,4)+n.xi_5*Math.pow(e,5)+n.xi_6*Math.pow(e,6)+n.xi_7*Math.pow(e,7)+n.xi_8*Math.pow(e,8)+n.xi_9*Math.pow(e,9)+n.xi_10*Math.pow(e,10)+n.xi_11*Math.pow(e,11)}function L(e,n){const _=Math.log(e),f=S(_,n);let a=Math.exp(G(_,n));return a=a<1?1:a,{sd:f,xi:a}}const C=2220446049250313e-31,H=5.656854249492381,J=.3989422804014327,V=[2.2352520354606837,161.02823106855587,1067.6894854603709,18154.98125334356,.06568233791820745],W=[47.202581904688245,976.0985517377767,10260.932208618979,45507.78933502673],X=[.39894151208813466,8.883149794388377,93.50665613217785,597.2702763948002,2494.5375852903726,6848.190450536283,11602.65143764735,9842.714838383978,10765576773720192e-24],Y=[22.266688044328117,235.387901782625,1519.3775994075547,6485.558298266761,18615.571640885097,34900.95272114598,38912.00328609327,19685.429676859992],x=[.215898534057957,.12740116116024736,.022235277870649807,.0014216191932278934,29112874951168793e-21,.023073441764940174],K=[1.284260096144911,.4682382124808651,.06598813786892856,.0037823963320275824,7297515550839662e-20],A=e=>e?-1/0:0,O=e=>e?0:1,R=(e,n)=>e?A(n):O(n),E=(e,n)=>e?O(n):A(n);function z(e,n,_,f,a){var r;if(isNaN(e)||isNaN(n)||isNaN(_))return NaN;if(!Number.isFinite(e)&&n===e)return NaN;if(_<=0)return _<0?NaN:e<n?R(f,a):E(f,a);if(r=(e-n)/_,!Number.isFinite(r))return e<n?R(f,a):E(f,a);e=r;var t=U(e,f?0:1,a);return f?t.cum:t.ccum}function w(e,n){return e===0?n?Number.NEGATIVE_INFINITY:0:n?0:1}function U(e,n,_){var f,a,r,t,s,i,o,c,l,u,p,g,m=V,N=W,h=X,v=Y,d=x,M=K;if(isNaN(e))return{cum:NaN,ccum:NaN};if(o=C*.5,p=n!==1,g=n!==0,l=Math.abs(e),l<=.67448975){if(l>o)for(c=e*e,t=m[4]*c,r=c,u=0;u<3;++u)t=(t+m[u])*c,r=(r+N[u])*c;else t=r=0;s=e*(t+m[3])/(r+N[3]),p&&(f=.5+s),g&&(a=.5-s),_&&(p&&(f=Math.log(f)),g&&(a=Math.log(a)))}else if(l<=H){for(t=h[8]*l,r=l,u=0;u<7;++u)t=(t+h[u])*l,r=(r+v[u])*l;s=(t+h[7])/(r+v[7]),c=Math.trunc(l*16)/16,i=(l-c)*(l+c),_?(f=-c*c*.5+-i*.5+Math.log(s),(p&&e>0||g&&e<=0)&&(a=Math.log1p(-Math.exp(-c*c*.5)*Math.exp(-i*.5)*s))):(f=Math.exp(-c*c*.5)*Math.exp(-i*.5)*s,a=1-f),e>0&&(s=f,p&&(f=a),a=s)}else if(_&&l<1e170||p&&-37.5193<e&&e<8.2924||g&&e<37.5193){for(c=1/(e*e),t=d[5]*c,r=c,u=0;u<4;++u)t=(t+d[u])*c,r=(r+M[u])*c;s=c*(t+d[4])/(r+M[4]),s=(J-s)/l,c=Math.trunc(e*16)/16,i=(e-c)*(e+c),_?(f=-c*c*.5+-i*.5+Math.log(s),(p&&e>0||g&&e<=0)&&(a=Math.log1p(-Math.exp(-c*c*.5)*Math.exp(-i*.5)*s))):(f=Math.exp(-c*c*.5)*Math.exp(-i*.5)*s,a=1-f),e>0&&(s=f,p&&(f=a),a=s)}else e>0?(f=w(1,_),a=w(0,_)):(f=w(0,_),a=w(1,_));return{cum:f,ccum:a}}})();