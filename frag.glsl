#version 430
#define debugmov 1
#define shadertoy 0 //noexport
#define doAA 0
#define partialrender 0
layout (location=0) uniform vec4 fpar[2];
layout (location=2) uniform vec4 debug[2]; //noexport
layout (location=4) uniform sampler2D tex;
#define MAT_GROUND 4
#define MAT_SCREEN_OUTER 5
#define MAT_WALL 6
#define MAT_PODIUM 7
#define MAT_OHP 8
#define MAT_OHP_LITE 9
#define MAT_BLU 10
int i;
vec3 gHitPosition = vec3(0);

mat2 rot2(float a){float s=sin(a),c=cos(a);return mat2(c,s,-s,c);}

// better box that works for subtraction //noexport
// https://iquilezles.org/articles/distfunctions/
float box(vec3 p, vec3 b)
{
	vec3 q = abs(p) - b;
	return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float ellipsoid(vec3 p, vec3 r)
{
  float k0 = length(p/r);
  float k1 = length(p/(r*r));
  return k0*(k0-1.0)/k1;
}

vec2 m(vec2 b, vec2 a){return a.x<b.x?a:b;}

vec2 wall(vec3 p)
{
	p.y -= 50.;
	return vec2(max(max(dot(p,vec3(0.,-1.,0.)), -box(p,vec3(60.,10.,68.))),
		-box(p,vec3(50.,30.,58.))
	), MAT_WALL);
}

// TODO: unused
vec2 screen(vec3 p)
{
	p.y -= 100.;
	float x = mod(p.x, 5.) - 2.5;
	return vec2(min(
		dot(p,vec3(0.,-1.,0.)), // wall
		length(max(abs(vec2(x, p.y)) - 1., 0.)) // bars
	), MAT_WALL);
}

// TODO: keep or not?
vec2 stage(vec3 p)
{
	p.y -= 68.;
	p.z += 1.;
	vec2 r = vec2(9e9, MAT_GROUND);
#define st r=m(r,vec2(length(max(abs(p) - vec3(8.,4.9,.15), 0.)), MAT_PODIUM))
	st;
	p.y += 10.;
	st;
	return r;
}

vec2 ohp(vec3 p)
{
	p.z += 5.;
	vec2 r = vec2(length(max(abs(p)-vec3(3.,3.,2.),0.))-.05, MAT_OHP);
	vec3 q=p;q.z-=.14;
	r=m(r,vec2(length(max(abs(q)-vec3(2.2),0.)),MAT_OHP_LITE));
	q=p+vec3(-3.25,2.,3.);
	r=m(r,vec2(length(max(abs(q)-vec3(.2,.2,4.),0.)), MAT_OHP)); // vertical beam thing
	r=m(r,vec2(length(max(abs(q+vec3(.1,-.4,2.5))-vec3(.5,.8,.5),0.)), MAT_OHP)); // handle thing
	q+=vec3(1.4,-.8,3.1);
	q.zx *= rot2(1.);
	q.zy *= rot2(-.3);
	r=m(r,vec2(length(max(abs(q)-vec3(.2,.2,1.7),0.)), MAT_OHP)); // connect thing
	p.z += 7.;
	r=m(r,vec2(max(dot(p-vec3(0.,0.,.3),vec3(0.,0.,-1.)),length(p)-.6),MAT_BLU)); // sphere
	r=m(r,vec2(max(
		length(max(abs(p)-vec3(.7,1.,.3),0.)), // lens
		dot(p,normalize(vec3(0.,1.,-2.))) // with slanted cutoff
	), MAT_OHP));
	p.z+=.8;
	p.yz*=rot2(-.9);
	r=m(r,vec2(length(max(abs(p)-vec3(.6,.8,.07),0.)), MAT_OHP));
	return r;
}

float beam(vec3 p)
{
	p.y-=35.;
	p.z+=28.; // 26+2 because map() does +2
	float t = dot(p+vec3(0.,0.,2.3),-vec3(0.,.43,1.)); // top
	t=max(t,dot(p+vec3(8.1,0.,0.),-vec3(1.,.22,0.))); // right
	t=max(t,dot(p-vec3(8.1,0.,0.),vec3(1.,-.22,0.))); // left
	t=max(t,dot(p-vec3(0.,0.,13.7),vec3(0.,0.,1.))); // bottom
	return max(t,dot(p+vec3(0.,35.,-13.2),vec3(0.,-.1,-0.08))); // back
}

vec2 map(vec3 p)
{
	p.z += 2.;

	vec2 r = vec2(9e9, MAT_GROUND);
	r = m(r, wall(p));
	//r = m(r, stage(p));
	r = m(r, ohp(p));
	return m(r, vec2(dot(p,vec3(0.,0.,-1.)), MAT_GROUND));
}

vec3 norm(vec3 p, float dist_to_p)
{
	vec2 e=vec2(1.,-1.)*.0035;
	return normalize(e.xyy*map(p+e.xyy).x+e.yyx*map(p+e.yyx).x+e.yxy*map(p+e.yxy).x+e.xxx*map(p+e.xxx).x);
}

int dc=1,d2=1;float cone=0.,c2=0.;
// x=hit(flopineShade) y=dist_to_p z=dist_to_ro w=material(if hit)
vec4 march(vec3 ro, vec3 rd)
{
	float b,dist;
	vec4 r = vec4(0);
	for (i = 0; i < 500 && r.z < 350.; i++){
		gHitPosition = ro + rd * r.z;
		vec2 m = map(gHitPosition);
		if (dc>0) {
			b = beam(gHitPosition);
			if (b<.0001) {
				dc=0;cone=.2;
			} else if (b<m.x) m=vec2(b,MAT_PODIUM);
		}
		if (d2>0) {
			b = length(max(abs(gHitPosition.xy)-vec2(2.2),0.));
			if (b<.0001) {
				d2=0;c2=clamp(.3+(gHitPosition.z+10.)/20.,0.,.2);
			} else if (b<m.x) m=vec2(b,MAT_PODIUM);
		}
		dist = m.x;
		if (dist < .0001 && m.y!=MAT_PODIUM) {
			r.x = float(i)/float(200); // TODO: this can just be 1. if not using flopine shade
			r.y = dist;
			r.w = m.y;
			break;
		}
		r.z += dist;
	}
	return r;
}

// sourced from https://www.shadertoy.com/view/lsKcDD
float calcAO(vec3 pos, vec3 nor)
{
	float occ = 0.;
	float sca = 1.;
	for(int i=0; i<5; i++) {
		float h = 0.001 + 0.15*float(i)/4.0;
		float d = map( pos + h*nor ).x;
		occ += (h-d)*sca;
		sca *= 0.95;
	}
	return clamp(1.0 - 1.5*occ, 0.0, 1.0);
}

// https://iquilezles.org/articles/rmshadows/
float softshadow(vec3 ro, vec3 rd)
{
	float res = 1.0;
	float ph = 9e9;
	for(float dist = .01; dist < 40.; ) {
		float h = map(ro + rd*dist).x;
		if (h<.001) {
			return 0.;
		}
		float y = h*h/(2.*ph);
		float d = sqrt(h*h-y*y);
		res = min(res, 128.*d/max(0.,dist-y)); // lower this number for more feather shadows
		ph = h;
		dist += h;
	}
	return res;
}

vec3 getmat(vec4 r)
{
	vec3 p = gHitPosition.xyz;
	switch (int(r.w)) {
	case MAT_SCREEN_OUTER: return vec3(.1);
	//case MAT_WALL: return vec3(222.,188.,153.)/255.;
	case MAT_WALL:
		if (p.y>70.) {
			if (abs(p.x)<35. && abs(p.z+32.)<20.) {
				return vec3(.5);
			}
			return vec3(.02);
		}
		//if (p.z<-60. && mod(p.x+.5,10.)<1.)
			//return vec3(146.,100.,64.)/255.;
		if (mod(p.z,10.)>1.)
			return vec3(233.,203.,169.)/255.;
		return vec3(224.,154.,69.)/255.;
	case MAT_PODIUM:
		return vec3(1.,.07,.07);
	case MAT_OHP:
		//return vec3(255.,220.,175.)/255.;
		return vec3(.2);
	case MAT_OHP_LITE:
		return vec3(1.);
	case MAT_BLU:
		return vec3(.2,.2,.8);
	}
	if (p.y>50.)
		return vec3(.01);
	return vec3(.53,.23,.09);
}

float mb(vec2 mb)
{
    mb.x=1.-mb.x;
    if (mb.x<0.||mb.y<0.||mb.x>1.||mb.y>1.)return 0.;
    float a=mod(mb.x,.04)<.02?0.:1.,b=1.;
    if (mod(mb.y,.04)<.02)a=1.-a;
    if (mb.x > .5) {
        b=a,a=1.;
    }
    mb.x = abs(mb.x-.5)*2.;
    if (mb.y > .33) {
        float y = (mb.y-.33)/.66;
        if (mb.x < .715) {
            float x = mb.x/.715;
            if (x > .5 && (x-.5)*2. > 1.-y) {
                return b;
            } else if (x > 1.-y) {
                return a;
            }
        } else if ((mb.x-.715)/.285<y) {
            return b;
        }
    }
    if (mb.x < mb.y) {
        return b;
    }
    return 0.;
}

float e(vec2 e,float w,float h){
	e*=e;
	return e.x/w+e.y/h;
}
float e2(vec2 z,float w,float h,float x,float y){
	z.x=abs(z.x)+x;z.y+=y;
	return e(z,w,h);
}
vec3 o(vec2 h)
{
	vec3 lit=.8*vec3(250.,255.,211.)/255.;
	h.y+=.8;
	if (abs(h.x)<2.&&abs(h.y)<1.) {
		lit*=clamp(mb((h+vec2(1.8,.8))/vec2(3.6,1.6)),.02,1.);
	}
	h.xy*=rot2(.4);
	h.y-=3.2;
	float z = e2(h,.1,.2,-.2,1.);
	if (z<1.) {// eyes
		if (.8<z) {
			return vec3(0.); // eye trace
		}
		if (e2(h,.06,.1,-.08,1.05)<.03) {
			return vec3(1.); // eye inner white
		}
		if (e2(h,.12,.4,-.12,1.)<.1) {
			return vec3(0.); // eye pupil black
		}
		return lit; // white
	}
	if (e2(h,2.4,.8,0.,.2)<.35) {
		// nose
		if (e2(h,1.,2.,-.4,.2)<.01) {
			return vec3(0.); // nose hole
		}
		return vec3(242.,158.,2.)/300.;
	}
	z = e2(h,.5,1.1,-.05,.6);
	if (z<1.) {
		if (.9<z) {
			return vec3(0.); // head trace
		}
		return lit;
	}
	h.x=abs(h.x)-.12;
	h.y+=.7;
		// gewei
		vec2 g=h*rot2(.15);
		g.y+=.9;
		if (e(g,.12,1.)<.2) {
			return vec3(242.,158.,2.)/300.;
		}
	h*=rot2(.6);
	h.y+=.9;
	z=e(h,.12,1.2);
	if (.12<z && z<.2) {
		return vec3(0.); // ears
	}
	h.y-=.2;
	if (e(h,.005,.7)<.2) {
		return vec3(0.); // ears inner
	}
	return lit;
}

vec3 colorHit(vec4 result, vec3 rd, vec3 normal, vec3 mat)
{
	//return mat;
	//return 1.-vec3(result.x);
	if (result.w == MAT_OHP_LITE) {
		// it's 4x4
		vec2 h=gHitPosition.xy;
		return o(h);
	}
	if (result.w == MAT_WALL) {
		vec2 h=gHitPosition.xz;
		h.y+=32.;
		if (abs(h.x)<18.&&abs(h.y)<18.) {
			return o(h/8.5);
		}
	}

	// https://www.shadertoy.com/view/lsKcDD
	// key light
	vec3 lig = normalize(vec3(-.2,-.6,-.15));

	// TODO HERE
	return mat * 2. * clamp(dot(normal, lig), .0, 1.);

	vec3 hal = normalize(lig-rd);
	float dif = clamp(dot(normal, lig), .0, 1.) * clamp(softshadow(gHitPosition, lig),.5,1.);

	float spe = pow(clamp(dot(normal, hal), .0, 1. ),16.0)* dif *
	(0.04 + 0.96*pow(clamp(1.0+dot(hal,rd),.0,1.), 5.0));

	vec3 col = mat * 1.*dif;
	// TODO: maybe no spe
	//col += 12.0*spe*vec3(1.0,0.7,0.5);

	// ambient light
	float occ = calcAO(gHitPosition, normal);
	float amb = clamp(0.5+0.5*normal.y, 0., 1.);
	col += mat*amb*occ*vec3(.12);

	// fog
	float t = result.z;
	//col *= exp(-0.00007*t*t);
	//col *= exp(-0.000007*t*t*t);
	//col *= exp(-0.00007*t);
        //col *= exp(-0.0005*t*t*t);
	return col;
}

out vec4 c;
in vec2 v;
void main()
{
	vec2 normuv = (v + 1.) / 2;
#if partialrender == 1
	vec2 from = fpar[0].zw;
	if (!(normuv.x < from.x && from.x < normuv.x + .11 &&
		normuv.y < from.y && from.y < normuv.y + .11))
	{
		c = vec4(texture(tex,normuv).xyz, 1.);
		return;
	}
#endif

	vec3 ro = vec3(-13., 1., -11.);
	vec3 at = vec3(-5., -10., 0.);

#if debugmov //noexport
	ro = debug[0].xyz/20.; //noexport
	float vertAngle = debug[1].y/20.; //noexport
	float horzAngle = debug[1].x/20.; //noexport
	if (abs(vertAngle) < .001) { //noexport
		vertAngle = .001; //noexport
	} //noexport
	float xylen = sin(vertAngle); //noexport
	vertAngle = cos(vertAngle); //noexport
	at.x = ro.x + cos(horzAngle) * xylen; //noexport
	at.y = ro.y + sin(horzAngle) * xylen; //noexport
	at.z = ro.z + vertAngle; //noexport
#endif //noexport

        vec3	cf = normalize(at-ro),
		cl = normalize(cross(cf,vec3(0,0,-1)));
	mat3 rdbase = mat3(cl,normalize(cross(cl,cf)),cf);

	vec3 resultcol = vec3(0.);
#if doAA == 1
	for (int aaa = 0; aaa < 2; aaa++) {
		for (int aab = 0; aab < 2; aab++) {
#else
	int aaa = 0, aab = 0;
#endif
#if shadertoy == 1 //noexport
			vec2 o = v + vec2(float(aab),float(aab)) / 2. - 0.5; //noexport
			vec2 uv = (o-.5*iResolution.xy)/iResolution.y; //noexport
#else //noexport
			vec2 iResolution = fpar[0].xy;
			vec2 uv = v*(iResolution + vec2(float(aaa),float(aab))/4)/iResolution;
			uv.y /= iResolution.x/iResolution.y;
#endif //noexport
			vec3 rd = rdbase*normalize(vec3(uv,1)), col = vec3(0.);

			vec4 result = march(ro, rd);

			if (result.x > 0.) { // hit
				vec3 normal = norm(gHitPosition, result.y);
				col = colorHit(result, rd, normal, getmat(result)*.3)+cone+c2;
			}
			resultcol += col;
#if doAA == 1
		}
	}
	resultcol /= 4.;
#endif

	c = vec4(pow(resultcol, vec3(.4545)), 1.0); // pow for gamma correction because all the cool kids do it
}
