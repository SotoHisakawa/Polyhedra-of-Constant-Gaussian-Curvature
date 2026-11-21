/* Polygon boundaries and all measurements refer to the paper's original vertices.
   Rendering triangles are private, including for nonconvex faces. */
(function (global) {
  "use strict";

  const sub = (a, b) => a.map((x, i) => x - b[i]);
  const dot = (a, b) => a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  const cross = (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
  const unit = a => { const n = Math.hypot(...a); return a.map(x => x/n); };
  const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
  const orient = (a, b, c) => (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);

  function faceNormal(vertices, face) {
    const origin = vertices[face[0]], area = [0, 0, 0];
    for (let j = 1; j < face.length-1; j++) {
      const wedge = cross(sub(vertices[face[j]], origin), sub(vertices[face[j+1]], origin));
      for (let i=0; i<3; i++) area[i] += wedge[i];
    }
    return unit(area);
  }

  function triangulate(vertices, face, normal) {
    if (face.length === 3) return [face];
    const origin = vertices[face[0]], e = unit(sub(vertices[face[1]], origin)), f = cross(normal, e);
    const xy = face.map(id => { const v = sub(vertices[id],origin); return [dot(v,e), dot(v,f)]; });
    const eps = 1e-12 * Math.max(...xy.map(v => dot2(v,v)));
    const loop = face.map((_, i) => i), triangles = [];
    while (loop.length > 3) {
      let found = false;
      for (let j=0; j<loop.length; j++) {
        const i0=loop[(j+loop.length-1)%loop.length], i1=loop[j], i2=loop[(j+1)%loop.length];
        const a=xy[i0], b=xy[i1], c=xy[i2];
        if (orient(a,b,c)<=eps) continue;
        const blocked=loop.some(i => i!==i0 && i!==i1 && i!==i2 &&
          orient(a,b,xy[i])>=-eps && orient(b,c,xy[i])>=-eps && orient(c,a,xy[i])>=-eps);
        if (blocked) continue;
        triangles.push([face[i0],face[i1],face[i2]]);
        loop.splice(j,1); found=true; break;
      }
      if (!found) throw new Error("Could not triangulate face "+face.map(i=>i+1).join("–"));
    }
    triangles.push(loop.map(i=>face[i]));
    return triangles;
  }

  const dot2 = (a,b) => a[0]*b[0] + a[1]*b[1];
  function boundaryEdges(mesh) {
    const edges = new Map();
    mesh.faces.forEach((face,faceId) => face.forEach((v,i) => {
      const w=face[(i+1)%face.length], a=Math.min(v,w), b=Math.max(v,w), key=a+","+b;
      if (!edges.has(key)) edges.set(key,{a,b,faces:[]});
      edges.get(key).faces.push(faceId);
    }));
    return [...edges.values()].sort((a,b)=>a.a-b.a || a.b-b.b);
  }

  function cornerRadians(mesh, faceId, i) {
    const face=mesh.faces[faceId], verts=mesh.vertices, pos=verts[face[i]];
    const incoming=sub(pos,verts[face[(i+face.length-1)%face.length]]);
    const outgoing=sub(verts[face[(i+1)%face.length]],pos);
    const normal=faceNormal(verts,face);
    return Math.PI-Math.atan2(dot(normal,cross(incoming,outgoing)),dot(incoming,outgoing));
  }

  function nearestOnSegment(a,b,x,y) {
    const dx=b[0]-a[0], dy=b[1]-a[1], length=dx*dx+dy*dy;
    const t=length ? clamp(((x-a[0])*dx+(y-a[1])*dy)/length,0,1) : 0;
    const px=a[0]+t*dx, py=a[1]+t*dy;
    return {x:px,y:py,z:a[2]+t*(b[2]-a[2]),distance:Math.hypot(px-x,py-y)};
  }

  class MeshViewer {
    static edges(mesh) { return boundaryEdges(mesh).map(({a,b})=>[a,b]); }
    static edgeLength(mesh,a,b) { return Math.hypot(...sub(mesh.vertices[a],mesh.vertices[b])); }
    static cornerAngle(mesh,face,corner) { return cornerRadians(mesh,face,corner)*180/Math.PI; }

    constructor(canvas, options={}) {
      this.canvas=canvas;
      this.context=canvas.getContext("2d");
      if (!this.context) throw new Error("Canvas 2D is unavailable");
      this.onPick=options.onPick || (()=>{});
      this.mode="none"; this.showEdges=true; this.showVertices=true; this.revealChanges=false;
      this.pointers=new Map(); this.listeners=[]; this.frame=0;
      this.pickTolerance=0.004;
      canvas.style.touchAction="none";
      const listen=(target,event,fn,opts) => {
        target.addEventListener(event,fn,opts);
        this.listeners.push(()=>target.removeEventListener(event,fn,opts));
      };
      listen(canvas,"pointerdown",e=>this.down(e));
      listen(canvas,"pointermove",e=>this.move(e));
      listen(canvas,"pointerup",e=>this.up(e,true));
      listen(canvas,"pointercancel",e=>this.up(e,false));
      listen(canvas,"lostpointercapture",e=>this.up(e,false));
      listen(canvas,"wheel",e=>{
        e.preventDefault();
        this.zoom(Math.exp(-clamp(e.deltaY * (e.deltaMode===1?16:1),-160,160)*0.0018));
      },{passive:false});
      listen(canvas,"keydown",e=>{
        const dy={ArrowUp:0.12,ArrowDown:-0.12}[e.key];
        const dx={ArrowLeft:-0.12,ArrowRight:0.12}[e.key];
        if (dx || dy) this.orbit(dx||0,dy||0);
        else if (e.key==="+" || e.key==="=") this.zoom(1.15);
        else if (e.key==="-") this.zoom(1/1.15);
        else if (e.key==="Home") this.reset();
        else return;
        e.preventDefault();
      });
      if (typeof ResizeObserver!=="undefined") {
        this.observer=new ResizeObserver(()=>this.schedule());
        this.observer.observe(canvas);
      } else listen(global,"resize",()=>this.schedule());
      this.reset();
    }

    destroy() {
      this.observer?.disconnect(); this.listeners.forEach(remove=>remove());
      if (this.frame) cancelAnimationFrame(this.frame);
    }

    setModel(mesh, changes={}, options={}) {
      const bounds=[0,1,2].map(k=>{
        const c=mesh.vertices.map(v=>v[k]); return (Math.min(...c)+Math.max(...c))/2;
      });
      const centred=mesh.vertices.map(v=>sub(v,bounds));
      const radius=Math.max(...centred.map(v=>Math.hypot(...v))) || 1;
      this.points=centred.map(v=>v.map(x=>x/radius));
      this.model=mesh;
      this.normals=mesh.faces.map(f=>faceNormal(mesh.vertices,f));
      this.triangles=mesh.faces.map((face,i)=>triangulate(mesh.vertices,face,this.normals[i]));
      this.edges=boundaryEdges(mesh);
      this.incident=mesh.vertices.map(()=>[]);
      mesh.faces.forEach((f,i)=>f.forEach(v=>this.incident[v].push(i)));
      this.changeFaces=new Set(changes?.faces || []);
      this.changeVertices=new Set(changes?.vertices || []);
      this.selection=null;
      this.measurementLabel=null;
      if (!options.preserveView) this.reset();
      this.schedule();
    }

    reset() { this.yaw=-0.63; this.pitch=0.53; this.magnification=1; this.schedule(); }
    zoom(factor) { if (factor>0) this.magnification=clamp(this.magnification*factor,0.3,5); this.schedule(); }
    orbit(dx,dy) { this.yaw+=dx; this.pitch=clamp(this.pitch+dy,-1.565,1.565); this.schedule(); }
    setWireframe(on) { this.showEdges=!!on; this.schedule(); }
    setVertices(on) { this.showVertices=!!on; this.schedule(); }
    setRevealChanges(on) { this.revealChanges=!!on; this.schedule(); }
    setPickMode(mode) { this.mode=mode; }
    setSelection(selection) { this.selection=selection; this.measurementLabel=null; this.schedule(); }
    setMeasurementLabel(text) {
      this.measurementLabel=text==null ? null : String(text);
      this.schedule();
    }
    schedule() {
      if (!this.frame) this.frame=requestAnimationFrame(()=>{this.frame=0;this.render();});
    }

    down(e) {
      if (e.button!==0 && e.pointerType!=="touch") return;
      this.canvas.setPointerCapture(e.pointerId);
      this.pointers.set(e.pointerId,{x:e.clientX,y:e.clientY,startX:e.clientX,startY:e.clientY});
      if (this.pointers.size===1) this.dragged=false;
      else this.dragged=true;
    }
    move(e) {
      const p=this.pointers.get(e.pointerId); if (!p) return;
      if (this.pointers.size===2) {
        const other=[...this.pointers.values()].find(v=>v!==p);
        const d0=Math.hypot(p.x-other.x,p.y-other.y);
        const d1=Math.hypot(e.clientX-other.x,e.clientY-other.y);
        if (d0>1 && d1>1) this.zoom(d1/d0);
      } else if (this.pointers.size===1) {
        if (Math.hypot(e.clientX-p.startX,e.clientY-p.startY)>4) this.dragged=true;
        if (this.dragged) this.orbit((e.clientX-p.x)*0.009,(e.clientY-p.y)*0.009);
      }
      p.x=e.clientX;p.y=e.clientY;
    }
    up(e,pick) {
      if (!this.pointers.delete(e.pointerId)) return;
      if (pick && !this.dragged && this.mode!=="none") this.selectAt(e.clientX,e.clientY);
    }

    prepareFrame() {
      const rect=this.canvas.getBoundingClientRect();
      if (!rect.width || !rect.height) return false;
      const dpi=Math.min(global.devicePixelRatio||1,1.6,1400/Math.max(rect.width,rect.height));
      const width=Math.max(1,Math.round(rect.width*dpi)), height=Math.max(1,Math.round(rect.height*dpi));
      this.dpi=width/rect.width;
      if (!this.image || width!==this.width || height!==this.height) {
        this.canvas.width=width; this.canvas.height=height;
        this.width=width; this.height=height;
        this.image=this.context.createImageData(width,height);
        this.pixels=new Uint32Array(this.image.data.buffer);
        this.depth=new Float32Array(width*height);
        this.faceAt=new Int32Array(width*height);
      }
      this.scale=Math.min(width,height)*0.44*this.magnification;
      this.cy=Math.cos(this.yaw); this.sy=Math.sin(this.yaw);
      this.cp=Math.cos(this.pitch); this.sp=Math.sin(this.pitch);
      this.pixels.fill(0);this.depth.fill(-Infinity);this.faceAt.fill(-1);
      return true;
    }
    rotate(v) {
      const q=this.sy*v[0]+this.cy*v[1];
      return [this.cy*v[0]-this.sy*v[1],this.cp*v[2]-this.sp*q,this.sp*v[2]+this.cp*q];
    }
    project(v) {
      const [x,y,z]=this.rotate(v);
      return [this.width/2+this.scale*x,this.height/2-this.scale*y,z];
    }

    render() {
      if (!this.model || !this.prepareFrame()) return;
      const projected=this.points.map(p=>this.project(p));
      this.projected=projected;
      // An initial base stage without highlighted features is shown in full.
      this.revealing=this.revealChanges && (this.changeFaces.size>0 || this.changeVertices.size>0);
      this.triangles.forEach((triangles,face)=>{
        const changed=this.changeFaces.has(face);
        if (this.revealing && !changed) return;
        const normal=this.rotate(this.normals[face]), sign=normal[2]>=0?1:-1;
        const shade=0.58+0.34*Math.max(0,sign*(-0.37*normal[0]+0.65*normal[1]+0.66*normal[2]))+0.08*Math.abs(normal[2]);
        const base=changed?[245,178,111]:[172,212,204];
        const col=base.map(v=>Math.round(v*shade));
        const packed=(255<<24)|(col[2]<<16)|(col[1]<<8)|col[0];
        for (const [i,j,k] of triangles) this.fillTriangle(projected[i],projected[j],projected[k],face,packed);
      });
      if (this.revealing && this.showEdges) for (const {a,b,faces} of this.edges) {
        if (!faces.some(id=>this.changeFaces.has(id))) {
          this.drawSegment(projected[a],projected[b],faces,[65,90,91],0.62*this.dpi,0.31);
        }
      }
      if (this.showEdges) for (const {a,b,faces} of this.edges) {
        if (this.revealing && !faces.some(id=>this.changeFaces.has(id))) continue;
        this.drawSegment(projected[a],projected[b],faces,[44,75,77],0.70*this.dpi);
      }
      if (this.showVertices || this.changeVertices.size) this.projected.forEach((p,i)=>{
        const marked=this.changeVertices.has(i);
        if (this.revealing && !marked && !this.incident[i].some(id=>this.changeFaces.has(id))) return;
        if (this.showVertices || marked) this.drawSegment(p,p,this.incident[i],marked?[174,75,32]:[51,82,81],(marked?3.0:1.65)*this.dpi);
      });
      this.context.putImageData(this.image,0,0);
      if (this.selection) this.drawSelection();
    }

    fillTriangle(a,b,c,face,color) {
      const area=orient(a,b,c);
      if (Math.abs(area)<1e-8) return;
      const xmin=Math.max(0,Math.floor(Math.min(a[0],b[0],c[0])));
      const xmax=Math.min(this.width-1,Math.ceil(Math.max(a[0],b[0],c[0])));
      const ymin=Math.max(0,Math.floor(Math.min(a[1],b[1],c[1])));
      const ymax=Math.min(this.height-1,Math.ceil(Math.max(a[1],b[1],c[1])));
      const px=xmin+0.5, py=ymin+0.5;
      const sx=(b[1]-c[1])/area, sy=(c[0]-b[0])/area;
      const tx=(c[1]-a[1])/area, ty=(a[0]-c[0])/area;
      let s0=((px-c[0])*(b[1]-c[1])+(py-c[1])*(c[0]-b[0]))/area;
      let t0=((px-c[0])*(c[1]-a[1])+(py-c[1])*(a[0]-c[0]))/area;
      const dzA=a[2]-c[2], dzB=b[2]-c[2];
      for (let y=ymin; y<=ymax; y++,s0+=sy,t0+=ty) {
        let s=s0,t=t0,index=y*this.width+xmin;
        for (let x=xmin; x<=xmax; x++,s+=sx,t+=tx,index++) {
          if (s< -1e-9 || t< -1e-9 || s+t>1+1e-9) continue;
          const z=c[2]+s*dzA+t*dzB;
          if (z>this.depth[index]) {
            this.depth[index]=z;this.faceAt[index]=face;this.pixels[index]=color;
          }
        }
      }
    }

    drawSegment(a,b,incident,color,radius,opacity=1) {
      const dx=b[0]-a[0],dy=b[1]-a[1],d2=dx*dx+dy*dy;
      const steps=Math.max(1,Math.ceil(Math.max(Math.abs(dx),Math.abs(dy))));
      const reach=Math.ceil(radius+0.5);
      // Traverse along the major axis so work is proportional to edge length.
      const horiz=Math.abs(dx)>=Math.abs(dy);
      const major=horiz?0:1, minor=1-major;
      const start=Math.max(0,Math.floor(Math.min(a[major],b[major]))-reach);
      const end=Math.min((horiz?this.width:this.height)-1,Math.ceil(Math.max(a[major],b[major]))+reach);
      for (let u=start;u<=end;u++) {
        const param=a[major]===b[major]?0:clamp((u+0.5-a[major])/(b[major]-a[major]),0,1);
        const mid=a[minor]+param*(b[minor]-a[minor]);
        const span=steps===1?reach:Math.ceil(1.42*(radius+0.5));
        const lo=Math.max(0,Math.floor(mid)-span), hi=Math.min((horiz?this.height:this.width)-1,Math.ceil(mid)+span);
        for (let v=lo;v<=hi;v++) {
          const x=horiz?u:v,y=horiz?v:u;
          const t=d2?clamp(((x+0.5-a[0])*dx+(y+0.5-a[1])*dy)/d2,0,1):0;
          const ex=x+0.5-a[0]-t*dx,ey=y+0.5-a[1]-t*dy;
          const cover=opacity*clamp(radius+0.5-Math.sqrt(ex*ex+ey*ey),0,1);
          if (!cover) continue;
          const i=y*this.width+x,z=a[2]+t*(b[2]-a[2]);
          if (z < this.depth[i]-this.pickTolerance && !incident.includes(this.faceAt[i])) continue;
          const p=4*i, data=this.image.data, alpha=data[p+3]/255, out=cover+alpha*(1-cover);
          for (let j=0;j<3;j++) data[p+j]=(color[j]*cover+data[p+j]*alpha*(1-cover))/out;
          data[p+3]=255*out;
        }
      }
    }

    selectAt(clientX,clientY) {
      if (this.frame || !this.projected) this.render();
      if (!this.projected) return;
      const r=this.canvas.getBoundingClientRect(),x=(clientX-r.left)*this.dpi,y=(clientY-r.top)*this.dpi;
      const at=(x,y)=>clamp(Math.floor(y),0,this.height-1)*this.width+clamp(Math.floor(x),0,this.width-1);
      let selected=null;
      if (this.mode==="edge") {
        const tolerance=11*this.dpi;
        let best=tolerance, bestPriority=-1;
        for (const {a,b,faces} of this.edges) {
          const p=nearestOnSegment(this.projected[a],this.projected[b],x,y), i=at(p.x,p.y);
          if (p.distance>tolerance || (p.z<this.depth[i]-this.pickTolerance && !faces.includes(this.faceAt[i]))) continue;
          const priority=Number(this.revealing && faces.some(id=>this.changeFaces.has(id)));
          if (priority<bestPriority || (priority===bestPriority && p.distance>best)) continue;
          best=p.distance;bestPriority=priority;selected={kind:"edge",a,b};
        }
      } else if (this.mode==="angle") {
        const faceId=this.faceAt[at(x,y)];
        if (faceId>=0) {
          const face=this.model.faces[faceId];
          let dist=Infinity,corner=0;
          face.forEach((id,j)=>{
            const p=this.projected[id],d=(p[0]-x)**2+(p[1]-y)**2;
            if(d<dist){dist=d;corner=j;}
          });
          selected={kind:"angle",face:faceId,corner};
        }
      }
      if (selected) { this.setSelection(selected);this.onPick(selected); }
    }

    drawSelection() {
      const ctx=this.context, sel=this.selection;
      const stroke=(points, colour, lineWidth)=>{
        ctx.strokeStyle=colour;ctx.lineWidth=lineWidth*this.dpi;ctx.beginPath();
        points.forEach(([x,y],i)=>i?ctx.lineTo(x,y):ctx.moveTo(x,y));ctx.stroke();
      };
      let path,arc=[];
      if (sel.kind==="edge" && this.projected[sel.a] && this.projected[sel.b]) {
        path=[this.projected[sel.a],this.projected[sel.b]];
      } else if (sel.kind==="angle") {
        const face=this.model.faces[sel.face];
        if (!face || !Number.isInteger(sel.corner) || sel.corner<0 || sel.corner>=face.length) return;
        const j=sel.corner, prev=face[(j+face.length-1)%face.length], cur=face[j], next=face[(j+1)%face.length];
        path=[this.projected[prev],this.projected[cur],this.projected[next]];
        const centre=this.points[cur],toPrev=sub(this.points[prev],centre),toNext=sub(this.points[next],centre);
        const radius=Math.min(Math.hypot(...toPrev),Math.hypot(...toNext))*0.27;
        const u=unit(toPrev),v=cross(u,this.normals[sel.face]);
        const theta=cornerRadians(this.model,sel.face,j);
        for(let k=0;k<=32;k++){
          const t=theta*k/32;
          arc.push(this.project(centre.map((x,i)=>x+radius*(Math.cos(t)*u[i]+Math.sin(t)*v[i]))));
        }
      } else return;
      ctx.save();ctx.lineCap="round";ctx.lineJoin="round";
      stroke(path,"rgba(255,255,255,.85)",4);stroke(path,"#ac4926",2);
      if(arc.length)stroke(arc,"#ac4926",2);
      for (const [x,y] of path) {
        ctx.beginPath();ctx.arc(x,y,3.2*this.dpi,0,2*Math.PI);
        ctx.fillStyle="#fff7e9";ctx.fill();ctx.lineWidth=1.5*this.dpi;ctx.strokeStyle="#a34122";ctx.stroke();
      }
      if (this.measurementLabel) {
        const anchor=sel.kind==="angle" ? path[1] : [
          (path[0][0]+path[1][0])/2, (path[0][1]+path[1][1])/2
        ];
        this.drawCallout(anchor,this.measurementLabel);
      }
      ctx.restore();
    }

    drawCallout(anchor,text) {
      const ctx=this.context, unit=this.dpi, margin=5*unit;
      const h=Math.min(24*unit,this.height-2*margin);
      ctx.font="600 "+(12*unit)+"px system-ui, sans-serif";
      const w=Math.min(ctx.measureText(text).width+16*unit,this.width-2*margin);
      if (w<=0 || h<=0) return;
      const x=clamp(anchor[0],margin,this.width-margin);
      const y=clamp(anchor[1],margin,this.height-margin);
      // Prefer a nearby patch of background over one which masks a face.
      let chosen=null, best=Infinity;
      for (const gap of [12*unit,29*unit]) {
        for (const left of [x+gap,x-gap-w]) for (const top of [y-gap-h,y+gap]) {
          const l=clamp(left,margin,this.width-margin-w), t=clamp(top,margin,this.height-margin-h);
          let covered=0;
          for (const sx of [.1,.35,.65,.9]) for (const sy of [.15,.5,.85]) {
            const row=Math.floor(t+sy*h),col=Math.floor(l+sx*w);
            if (this.faceAt[row*this.width+col]>=0) covered++;
          }
          const px=clamp(x,l,l+w),py=clamp(y,t,t+h);
          const distance=Math.hypot(px-x,py-y);
          const score=covered+distance/(20*unit)+(distance<4*unit?20:0);
          if (score<best) { best=score;chosen={l,t,px,py}; }
        }
      }
      const {l,t,px,py}=chosen;
      ctx.beginPath();ctx.moveTo(x,y);ctx.lineTo(px,py);
      ctx.lineWidth=unit;ctx.strokeStyle="#924322";ctx.stroke();
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(l,t,w,h,Math.min(4*unit,h/2,w/2));
      else ctx.rect(l,t,w,h);
      ctx.fillStyle="#fffaf2";ctx.fill();ctx.strokeStyle="#b27d61";ctx.stroke();
      ctx.fillStyle="#382d27";ctx.textBaseline="middle";ctx.textAlign="center";
      ctx.fillText(text,l+w/2,t+h/2,Math.max(1,w-9*unit));
    }
  }
  global.MeshViewer=MeshViewer;
})(typeof window === "undefined" ? globalThis : window);
