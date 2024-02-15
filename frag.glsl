#version 430
#define debugmov 1
#define shadertoy 0 //noexport
#define doAA 0
#define partialrender 0
layout (location=0) uniform vec4 fpar[2];
layout (location=2) uniform vec4 debug[2]; //noexport
layout (location=4) uniform sampler2D tex;
#define MAT_KEY_BLACK 0
#define MAT_KEY_WHITE 1
#define MAT_BLACK_NOISE 2
#define MAT_BLACK_NOISE_LOGO 3
#define MAT_GROUND 4
#define MAT_SCREEN_OUTER 5
#define MAT_WALL 6
int i;
vec3 gHitPosition = vec3(0);

float rand(vec2 p){return fract(sin(dot(p.xy,vec2(12.9898,78.233)))*43758.5453);}
mat2 rot2(float a){float s=sin(a),c=cos(a);return mat2(c,s,-s,c);}

// better box that works for subtraction //noexport
// https://iquilezles.org/articles/distfunctions/
float box(vec3 p, vec3 b)
{
	vec3 q = abs(p) - b;
	return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

vec2 m(vec2 b, vec2 a){return a.x<b.x?a:b;}

vec2 wall(vec3 p)
{
	p.y -= 50.;
	return vec2(max(max(dot(p,vec3(0.,-1.,0.)), -box(p,vec3(60.,10.,68.))),
		-box(p,vec3(50.,30.,58.))
	), MAT_WALL);
}

vec2 screen(vec3 p)
{
	p.y -= 100.;
	float x = mod(p.x, 5.) - 2.5;
	return vec2(min(
		dot(p,vec3(0.,-1.,0.)), // wall
		length(max(abs(vec2(x, p.y)) - 1., 0.)) // bars
	), MAT_WALL);
}

vec2 map(vec3 p)
{
	float ground = dot(p,vec3(0.,0.,-1.));
	p.z += 2.;

	vec2 r = vec2(9e9, MAT_GROUND);
	r = m(r, vec2(length(max(abs(p.xyz-vec3(8.,0.,0.)) - vec3(.4,2.,1.), 0.)), MAT_SCREEN_OUTER));
	r = m(r, vec2(length(max(abs(p.yz-vec2(8.1,2.4)) - vec2(.4,2.), 0.)), MAT_BLACK_NOISE));
	r = m(r, wall(p));
	if (ground < r.x) return vec2(ground, MAT_GROUND);
	return r;
}

vec3 norm(vec3 p, float dist_to_p)
{
	vec2 e=vec2(1.,-1.)*.0035;
	return normalize(e.xyy*map(p+e.xyy).x+e.yyx*map(p+e.yyx).x+e.yxy*map(p+e.yxy).x+e.xxx*map(p+e.xxx).x);
}

// x=hit(flopineShade) y=dist_to_p z=dist_to_ro w=material(if hit)
vec4 march(vec3 ro, vec3 rd)
{
	vec4 r = vec4(0);
	for (i = 0; i < 200 && r.z < 350.; i++){
		gHitPosition = ro + rd * r.z;
		vec2 m = map(gHitPosition);
		float dist = m.x;
		if (dist < .0001) {
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

// w component is amount of reflection mix //noexport
vec4 getmat(vec4 r)
{
	vec3 p = gHitPosition.xyz;
	switch (int(r.w)) {
	case MAT_KEY_BLACK: return vec4(.007,.007,.007,.4);
	case MAT_KEY_WHITE: return vec4(vec3(218.,216.,227.)/255., .6);
	case MAT_BLACK_NOISE:
	case MAT_BLACK_NOISE_LOGO: return vec4(vec3(.05+.05*rand(mod(vec2(r.z,r.y),10))), 0.);
	case MAT_GROUND:
		if (p.y>50.)
			return vec4(.01);
		return vec4(.53,.23,.09, 0.);
	case MAT_SCREEN_OUTER: return vec4(vec3(.1),.0);
	//case MAT_WALL: return vec4(222.,188.,153.,.0)/255.;
	case MAT_WALL:
		if (p.y>70.) {
			if (abs(p.x)<35. && abs(p.z+32.)<20.) {
				return vec4(vec3(.5),0.);
			}
			return vec4(vec3(.02),0.);
		}
		//if (p.z<-60. && mod(p.x+.5,10.)<1.)
			//return vec4(146.,100.,64.,.0)/255.;
		if (mod(p.z,10.)>1.)
			return vec4(233.,203.,169.,.0)/255.;
		return vec4(224.,154.,69.,0.)/255.;
	}
	return vec4(0., 1., 0., 3.);
}

vec3 colorHit(vec4 result, vec3 rd, vec3 normal, vec3 mat)
{
	//return mat;
	//return 1.-vec3(result.x);
	if (result.w == MAT_BLACK_NOISE_LOGO) {
		vec2 mb = gHitPosition.xy;
		// x-4 is mid, height is 4 so x is -8 to 0
		if (mb.x > -8. && mb.x < 0. && mb.y < -18.7 && mb.y > -22.7) {
			mb -= vec2(-8., -22.7);
			mb /= vec2(8., 4.);
			float a = .2, b = .8;
			if (mb.x > .5) {
				b = .2;a = .8;
			}
			mb.x = abs(mb.x-.5)*2.;
			if (mb.y > .33) {
				float y = (mb.y-.33)/.66;
				if (mb.x < .715) {
					float x = mb.x/.715;
					if (x > .5 && (x-.5)*2. > 1.-y) {
						return vec3(1.)*b;
					} else if (x > 1.-y) {
						return vec3(1.)*a;
					}
				} else if ((mb.x-.715)/.285<y) {
					return vec3(1.)*b;
				}
			}
			if (mb.x < mb.y) {
				return vec3(1.)*b;
			}
		}
	}

	// https://www.shadertoy.com/view/lsKcDD
	// key light
	vec3 lig = normalize(vec3(-.2,-.6,-.15));
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
	col += mat*amb*occ*vec3(.1);

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
				vec4 mat = getmat(result) * .3;
				// reflexxions
				if (mat.w > .0001) {
					vec3 gg = gHitPosition;
					rd = reflect(rd, normal);
					gHitPosition += .001 * rd;
					vec4 nr = march(gHitPosition, rd);
					if (result.x > 0.) {
						vec3 nn = norm(gHitPosition, result.y);
						vec3 m = getmat(nr).xyz;
						vec3 nc = colorHit(nr, rd, nn, m);
						mat.xyz = mix(mat.xyz, nc * mat.w, mat.w);
					}
					gHitPosition = gg;
				}
				col = colorHit(result, rd, normal, mat.xyz);
			}
			resultcol += col;
#if doAA == 1
		}
	}
	resultcol /= 4.;
#endif

	c = vec4(pow(resultcol, vec3(.4545)), 1.0); // pow for gamma correction because all the cool kids do it
}
