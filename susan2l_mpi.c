/* }}} */
/* {{{ Readme First */

/**********************************************************************
 SUSAN Version 2l
 SUSAN = Smallest Univalue Segment Assimilating Nucleus

 Compile with:
 gcc -O4 -o susan susan2l.c -lm

 See following section for different machine information. Please
 report any bugs (and fixes). There are a few optional changes that
 can be made in the "defines" section which follows shortly.

 Usage: type "susan" to get usage. Only PGM format files can be input
 and output. Utilities such as the netpbm package and XV can be used
 to convert to and from other formats. Any size of image can be
 processed.

 This code is written using an emacs folding mode, making moving
 around the different sections very easy. This is why there are
 various marks within comments and why comments are indented.


 SUSAN QUICK:

 This version of the SUSAN corner finder does not do all the
 false-corner suppression and thus is faster and produced some false
 positives, particularly on strong edges. However, because there are
 less stages involving thresholds etc., the corners that are
 correctly reported are usually more stable than those reported with
 the full algorithm. Thus I recommend at least TRYING this algorithm
 for applications where stability is important, e.g., tracking.

 THRESHOLDS:

 There are two thresholds which can be set at run-time. These are the
 brightness threshold (t) and the distance threshold (d).

 SPATIAL CONTROL: d

 In SUSAN smoothing d controls the size of the Gaussian mask; its
 default is 4.0. Increasing d gives more smoothing. In edge finding,
 a fixed flat mask is used, either 37 pixels arranged in a "circle"
 (default), or a 3 by 3 mask which gives finer detail. In corner
 finding, only the larger 37 pixel mask is used; d is not
 variable. In smoothing, the flat 3 by 3 mask can be used instead of
 a larger Gaussian mask; this gives low smoothing and fast operation.

 BRIGHTNESS CONTROL: t

 In all three algorithms, t can be varied (default=20); this is the
 main threshold to be varied. It determines the maximum difference in
 greylevels between two pixels which allows them to be considered
 part of the same "region" in the image. Thus it can be reduced to
 give more edges or corners, i.e. to be more sensitive, and vice
 versa. In smoothing, reducing t gives less smoothing, and vice
 versa. Set t=10 for the test image available from the SUSAN web
 page.

 ITERATIONS:

 With SUSAN smoothing, more smoothing can also be obtained by
 iterating the algorithm several times. This has a different effect
 from varying d or t.

 FIXED MASKS:

 37 pixel mask:    ooo       3 by 3 mask:  ooo
 ooooo                    ooo
 ooooooo                   ooo
 ooooooo
 ooooooo
 ooooo
 ooo

 CORNER ATTRIBUTES dx, dy and I
 (Only read this if you are interested in the C implementation or in
 using corner attributes, e.g., for corner matching)

 Corners reported in the corner list have attributes associated with
 them as well as positions. This is useful, for example, when
 attempting to match corners from one image to another, as these
 attributes can often be fairly unchanged between images. The
 attributes are dx, dy and I. I is the value of image brightness at
 the position of the corner. In the case of susan_corners_quick, dx
 and dy are the first order derivatives (differentials) of the image
 brightness in the x and y directions respectively, at the position
 of the corner. In the case of normal susan corner finding, dx and dy
 are scaled versions of the position of the centre of gravity of the
 USAN with respect to the centre pixel (nucleus).

 BRIGHTNESS FUNCTION LUT IMPLEMENTATION:
 (Only read this if you are interested in the C implementation)

 The SUSAN brightness function is implemented as a LUT
 (Look-Up-Table) for speed. The resulting pointer-based code is a
 little hard to follow, so here is a brief explanation. In
 setup_brightness_lut() the LUT is setup. This mallocs enough space
 for *bp and then repositions the pointer to the centre of the
 malloced space. The SUSAN function e^-(x^6) or e^-(x^2) is
 calculated and converted to a uchar in the range 0-100, for all
 possible image brightness differences (including negative
 ones). Thus bp[23] is the output for a brightness difference of 23
 greylevels. In the SUSAN algorithms this LUT is used as follows:

 p=in + (i-3)*x_size + j - 1;
 p points to the first image pixel in the circular mask surrounding
 point (x,y).

 cp=bp + in[i*x_size+j];
 cp points to a position in the LUT corresponding to the brightness
 of the centre pixel (x,y).

 now for every pixel within the mask surrounding (x,y),
 n+=*(cp-*p++);
 the brightness difference function is found by moving the cp pointer
 down by an amount equal to the value of the pixel pointed to by p,
 thus subtracting the two brightness values and performing the
 exponential function. This value is added to n, the running USAN
 area.

 in SUSAN smoothing, the variable height mask is implemented by
 multiplying the above by the moving mask pointer, reset for each new
 centre pixel.
 tmp = *dpt++ * *(cp-brightness);

 \**********************************************************************/

/* }}} */
/* {{{ defines, includes and typedefs */

/* ********** Optional settings */

#ifndef PPC
typedef int TOTAL_TYPE; /* this is faster for "int" but should be "float" for large d masks */
#else
typedef float TOTAL_TYPE; /* for my PowerPC accelerator only */
#endif

/*#define FOPENB*//* uncomment if using djgpp gnu C for DOS or certain Win95 compilers */
#define SEVEN_SUPP           /* size for non-max corner suppression; SEVEN_SUPP or FIVE_SUPP */
#define MAX_CORNERS   15000  /* max corners per frame */

/* ********** Leave the rest - but you may need to remove one or both of sys/file.h and malloc.h lines */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/file.h>    /* may want to remove this line */
#include <malloc.h>      /* may want to remove this line */
#include <mpi.h>

#define  	MASTER        0          /* process ID for MASTER */
#define 	BEGIN         1          /* Message tag */
#define 	DONE          2          /* Message tag */

#define  	exit_error(IFB,IFC) { fprintf(stderr,IFB,IFC); exit(0); }
#define  	FTOI(a) ( (a) < 0 ? ((int)(a-0.5)) : ((int)(a+0.5)) )

typedef unsigned char uchar;

/* }}} */
/* {{{ get_image(filename,in,x_size,y_size) */

/* {{{ int getint(fp) derived from XV */

int getint(fd)
	FILE *fd; {
	int c, i;
	char dummy[10000];

	c = getc(fd);
	while (1) /* find next integer */
	{
		if (c == '#') /* if we're at a comment, read to end of line */
			fgets(dummy, 9000, fd);
		if (c == EOF)
			exit_error("Image %s not binary PGM.\n", "is");
		if (c >= '0' && c <= '9')
			break; /* found what we were looking for */
		c = getc(fd);
	}

	/* we're at the start of a number, continue until we hit a non-number */
	i = 0;
	while (1) {
		i = (i * 10) + (c - '0');
		c = getc(fd);
		if (c == EOF)
			return (i);
		if (c < '0' || c > '9')
			break;
	}

	return (i);
}

/* }}} */

/* reads whole image in "in" variable" */
void get_image(filename, in, x_size, y_size)
	char filename[200];unsigned char **in;int *x_size, *y_size; {
	FILE *fd;
	char header[100];
	int tmp;

#ifdef FOPENB
	if ((fd=fopen(filename,"rb")) == NULL)
#else
	if ((fd = fopen(filename, "r")) == NULL)
#endif
		exit_error("Can't input image %s.\n", filename);

	/* {{{ read header */

	header[0] = fgetc(fd);
	header[1] = fgetc(fd);
	if (!(header[0] == 'P' && header[1] == '5'))
		exit_error("Image %s does not have binary PGM header.\n", filename);

	*x_size = getint(fd);
	*y_size = getint(fd);
	tmp = getint(fd);

	/* }}} */

	*in = (uchar *) malloc(*x_size * *y_size);

	if (fread(*in, 1, *x_size * *y_size, fd) == 0)
		exit_error("Image %s is wrong size.\n", filename);

	fclose(fd);
}

/* }}} */
/* {{{ put_image(filename,in,x_size,y_size) */

put_image(filename, in, x_size, y_size)
	char filename[100], *in;int x_size, y_size; {
	FILE *fd;

#ifdef FOPENB
	if ((fd=fopen(filename,"wb")) == NULL)
#else
	if ((fd = fopen(filename, "w")) == NULL)
#endif
		exit_error("Can't output image%s.\n", filename);

	fprintf(fd, "P5\n");
	fprintf(fd, "%d %d\n", x_size, y_size);
	fprintf(fd, "255\n");

	if (fwrite(in, x_size * y_size, 1, fd) != 1)
		exit_error("Can't write image %s.\n", filename);

	fclose(fd);
}

/* }}} */
/* {{{ int_to_uchar(r,in,size) */

int_to_uchar(r, in, size)
	uchar *in;int *r, size; {
	int i, max_r = r[0], min_r = r[0];

	for (i = 0; i < size; i++) {
		if (r[i] > max_r)
			max_r = r[i];
		if (r[i] < min_r)
			min_r = r[i];
	}

	/*printf("min=%d max=%d\n",min_r,max_r);*/

	max_r -= min_r;

	//TODO
	for (i = 0; i < size; i++)
		in[i] = (uchar) ((int) ((int) (r[i] - min_r) * 255) / max_r);
}

/* }}} */
/* {{{ setup_brightness_lut(bp,thresh,form) */

void setup_brightness_lut(bp, thresh, form)
	uchar **bp;int thresh, form; {
	int k;
	float temp;

	*bp = (unsigned char *) malloc(516);
	*bp = *bp + 258;
	
	//TODO
	for (k = -256; k < 257; k++) {
		temp = ((float) k) / ((float) thresh);
		temp = temp * temp;
		if (form == 6)
			temp = temp * temp * temp;
		temp = 100.0 * exp(-temp);
		*(*bp + k) = (uchar) temp;
	}
}

/* }}} */

/* {{{ edges */

/* {{{ edge_draw(in,corner_list,drawing_mode) */

edge_draw(in, mid, x_size, y_size, drawing_mode)
	uchar *in, *mid;int x_size, y_size, drawing_mode; {
	int i;
	uchar *inp, *midp;
	
	//TODO
	if (drawing_mode == 0) {
		/* mark 3x3 white block around each edge point */
		midp = mid;
		for (i = 0; i < x_size * y_size; i++) {
			if (*midp < 8) {
				inp = in + (midp - mid) - x_size - 1;
				*inp++ = 255;
				*inp++ = 255;
				*inp = 255;
				inp += x_size - 2;
				*inp++ = 255;
				*inp++;
				*inp = 255;
				inp += x_size - 2;
				*inp++ = 255;
				*inp++ = 255;
				*inp = 255;
			}
			midp++;
		}
	}

	/* now mark 1 black pixel at each edge point */
	midp = mid;
	for (i = 0; i < x_size * y_size; i++) {
		if (*midp < 8)
			*(in + (midp - mid)) = 0;
		midp++;
	}
}

/* }}} */
/* {{{ susan_thin(r,mid,x_size,y_size) */

/* only one pass is needed as i,j are decremented if necessary to go
 back and do bits again */

susan_thin(r, mid, x_size, y_size)
	uchar *mid;int *r, x_size, y_size; {
	int l[9], centre, nlinks, npieces, b01, b12, b21, b10, p1, p2, p3, p4, b00,
			b02, b20, b22, m, n, a, b, x, y, i, j;
	uchar *mp;


	//TODO
	for (i = 4; i < y_size - 4; i++)
		for (j = 4; j < x_size - 4; j++)
			if (mid[i * x_size + j] < 8) {
				centre = r[i * x_size + j];
				/* {{{ count number of neighbours */

				mp = mid + (i - 1) * x_size + j - 1;

				n = (*mp < 8) + (*(mp + 1) < 8) + (*(mp + 2) < 8)
						+ (*(mp + x_size) < 8) + (*(mp + x_size + 2) < 8)
						+ (*(mp + x_size + x_size) < 8)
						+ (*(mp + x_size + x_size + 1) < 8)
						+ (*(mp + x_size + x_size + 2) < 8);

				/* }}} */
				/* {{{ n==0 no neighbours - remove point */

				if (n == 0)
					mid[i * x_size + j] = 100;

				/* }}} */
				/* {{{ n==1 - extend line if I can */

				/* extension is only allowed a few times - the value of mid is used to control this */

				if ((n == 1) && (mid[i * x_size + j] < 6)) {
					/* find maximum neighbour weighted in direction opposite the
					 neighbour already present. e.g.
					 have: O O O  weight r by 0 2 3
					 X X O              0 0 4
					 O O O              0 2 3     */

					l[0] = r[(i - 1) * x_size + j - 1];
					l[1] = r[(i - 1) * x_size + j];
					l[2] = r[(i - 1) * x_size + j + 1];
					l[3] = r[(i) * x_size + j - 1];
					l[4] = 0;
					l[5] = r[(i) * x_size + j + 1];
					l[6] = r[(i + 1) * x_size + j - 1];
					l[7] = r[(i + 1) * x_size + j];
					l[8] = r[(i + 1) * x_size + j + 1];

					if (mid[(i - 1) * x_size + j - 1] < 8) {
						l[0] = 0;
						l[1] = 0;
						l[3] = 0;
						l[2] *= 2;
						l[6] *= 2;
						l[5] *= 3;
						l[7] *= 3;
						l[8] *= 4;
					} else {
						if (mid[(i - 1) * x_size + j] < 8) {
							l[1] = 0;
							l[0] = 0;
							l[2] = 0;
							l[3] *= 2;
							l[5] *= 2;
							l[6] *= 3;
							l[8] *= 3;
							l[7] *= 4;
						} else {
							if (mid[(i - 1) * x_size + j + 1] < 8) {
								l[2] = 0;
								l[1] = 0;
								l[5] = 0;
								l[0] *= 2;
								l[8] *= 2;
								l[3] *= 3;
								l[7] *= 3;
								l[6] *= 4;
							} else {
								if (mid[(i) * x_size + j - 1] < 8) {
									l[3] = 0;
									l[0] = 0;
									l[6] = 0;
									l[1] *= 2;
									l[7] *= 2;
									l[2] *= 3;
									l[8] *= 3;
									l[5] *= 4;
								} else {
									if (mid[(i) * x_size + j + 1] < 8) {
										l[5] = 0;
										l[2] = 0;
										l[8] = 0;
										l[1] *= 2;
										l[7] *= 2;
										l[0] *= 3;
										l[6] *= 3;
										l[3] *= 4;
									} else {
										if (mid[(i + 1) * x_size + j - 1] < 8) {
											l[6] = 0;
											l[3] = 0;
											l[7] = 0;
											l[0] *= 2;
											l[8] *= 2;
											l[1] *= 3;
											l[5] *= 3;
											l[2] *= 4;
										} else {
											if (mid[(i + 1) * x_size + j] < 8) {
												l[7] = 0;
												l[6] = 0;
												l[8] = 0;
												l[3] *= 2;
												l[5] *= 2;
												l[0] *= 3;
												l[2] *= 3;
												l[1] *= 4;
											} else {
												if (mid[(i + 1) * x_size + j + 1]
														< 8) {
													l[8] = 0;
													l[5] = 0;
													l[7] = 0;
													l[6] *= 2;
													l[2] *= 2;
													l[1] *= 3;
													l[3] *= 3;
													l[0] *= 4;
												}
											}
										}
									}
								}
							}
						}
					}

					m = 0; /* find the highest point */
					for (y = 0; y < 3; y++)
						for (x = 0; x < 3; x++)
							if (l[y + y + y + x] > m) {
								m = l[y + y + y + x];
								a = y;
								b = x;
							}

					if (m > 0) {
						if (mid[i * x_size + j] < 4)
							mid[(i + a - 1) * x_size + j + b - 1] = 4;
						else
							mid[(i + a - 1) * x_size + j + b - 1] = mid[i
									* x_size + j] + 1;
						if ((a + a + b) < 3) /* need to jump back in image */
						{
							i += a - 1;
							j += b - 2;
							if (i < 4)
								i = 4;
							if (j < 4)
								j = 4;
						}
					}
				}

				/* }}} */
				/* {{{ n==2 */

				if (n == 2) {
					/* put in a bit here to straighten edges */
					b00 = mid[(i - 1) * x_size + j - 1] < 8; /* corners of 3x3 */
					b02 = mid[(i - 1) * x_size + j + 1] < 8;
					b20 = mid[(i + 1) * x_size + j - 1] < 8;
					b22 = mid[(i + 1) * x_size + j + 1] < 8;
					if (((b00 + b02 + b20 + b22) == 2)
							&& ((b00 | b22) & (b02 | b20))) { /* case: move a point back into line.
							 e.g. X O X  CAN  become X X X
							 O X O              O O O
							 O O O              O O O    */
						if (b00) {
							if (b02) {
								x = 0;
								y = -1;
							} else {
								x = -1;
								y = 0;
							}
						} else {
							if (b02) {
								x = 1;
								y = 0;
							} else {
								x = 0;
								y = 1;
							}
						}
						if (((float) r[(i + y) * x_size + j + x]
								/ (float) centre) > 0.7) {
							if (((x == 0)
									&& (mid[(i + (2 * y)) * x_size + j] > 7)
									&& (mid[(i + (2 * y)) * x_size + j - 1] > 7)
									&& (mid[(i + (2 * y)) * x_size + j + 1] > 7))
									|| ((y == 0)
											&& (mid[(i) * x_size + j + (2 * x)]
													> 7)
											&& (mid[(i + 1) * x_size + j
													+ (2 * x)] > 7)
											&& (mid[(i - 1) * x_size + j
													+ (2 * x)] > 7))) {
								mid[(i) * x_size + j] = 100;
								mid[(i + y) * x_size + j + x] = 3; /* no jumping needed */
							}
						}
					} else {
						b01 = mid[(i - 1) * x_size + j] < 8;
						b12 = mid[(i) * x_size + j + 1] < 8;
						b21 = mid[(i + 1) * x_size + j] < 8;
						b10 = mid[(i) * x_size + j - 1] < 8;
						/* {{{ right angle ends - not currently used */

#ifdef IGNORETHIS
						if ( (b00&b01)|(b00&b10)|(b02&b01)|(b02&b12)|(b20&b10)|(b20&b21)|(b22&b21)|(b22&b12) )
						{ /* case; right angle ends. clean up.
						 e.g.; X X O  CAN  become X X O
						 O X O              O O O
						 O O O              O O O        */
							if ( ((b01)&(mid[(i-2)*x_size+j-1]>7)&(mid[(i-2)*x_size+j]>7)&(mid[(i-2)*x_size+j+1]>7)&
											((b00&((2*r[(i-1)*x_size+j+1])>centre))|(b02&((2*r[(i-1)*x_size+j-1])>centre)))) |
									((b10)&(mid[(i-1)*x_size+j-2]>7)&(mid[(i)*x_size+j-2]>7)&(mid[(i+1)*x_size+j-2]>7)&
											((b00&((2*r[(i+1)*x_size+j-1])>centre))|(b20&((2*r[(i-1)*x_size+j-1])>centre)))) |
									((b12)&(mid[(i-1)*x_size+j+2]>7)&(mid[(i)*x_size+j+2]>7)&(mid[(i+1)*x_size+j+2]>7)&
											((b02&((2*r[(i+1)*x_size+j+1])>centre))|(b22&((2*r[(i-1)*x_size+j+1])>centre)))) |
									((b21)&(mid[(i+2)*x_size+j-1]>7)&(mid[(i+2)*x_size+j]>7)&(mid[(i+2)*x_size+j+1]>7)&
											((b20&((2*r[(i+1)*x_size+j+1])>centre))|(b22&((2*r[(i+1)*x_size+j-1])>centre)))) )
							{
								mid[(i)*x_size+j]=100;
								if (b10&b20) j-=2;
								if (b00|b01|b02) {i--; j-=2;}
							}
						}
#endif

						/* }}} */
						if (((b01 + b12 + b21 + b10) == 2)
								&& ((b10 | b12) & (b01 | b21))
								&& ((b01
										& ((mid[(i - 2) * x_size + j - 1] < 8)
												| (mid[(i - 2) * x_size + j + 1]
														< 8)))
										| (b10
												& ((mid[(i - 1) * x_size + j - 2]
														< 8)
														| (mid[(i + 1) * x_size
																+ j - 2] < 8)))
										| (b12
												& ((mid[(i - 1) * x_size + j + 2]
														< 8)
														| (mid[(i + 1) * x_size
																+ j + 2] < 8)))
										| (b21
												& ((mid[(i + 2) * x_size + j - 1]
														< 8)
														| (mid[(i + 2) * x_size
																+ j + 1] < 8))))) { /* case; clears odd right angles.
								 e.g.; O O O  becomes O O O
								 X X O          X O O
								 O X O          O X O     */
							mid[(i) * x_size + j] = 100;
							i--; /* jump back */
							j -= 2;
							if (i < 4)
								i = 4;
							if (j < 4)
								j = 4;
						}
					}
				}

				/* }}} */
				/* {{{ n>2 the thinning is done here without breaking connectivity */

				if (n > 2) {
					b01 = mid[(i - 1) * x_size + j] < 8;
					b12 = mid[(i) * x_size + j + 1] < 8;
					b21 = mid[(i + 1) * x_size + j] < 8;
					b10 = mid[(i) * x_size + j - 1] < 8;
					if ((b01 + b12 + b21 + b10) > 1) {
						b00 = mid[(i - 1) * x_size + j - 1] < 8;
						b02 = mid[(i - 1) * x_size + j + 1] < 8;
						b20 = mid[(i + 1) * x_size + j - 1] < 8;
						b22 = mid[(i + 1) * x_size + j + 1] < 8;
						p1 = b00 | b01;
						p2 = b02 | b12;
						p3 = b22 | b21;
						p4 = b20 | b10;

						if (((p1 + p2 + p3 + p4)
								- ((b01 & p2) + (b12 & p3) + (b21 & p4)
										+ (b10 & p1))) < 2) {
							mid[(i) * x_size + j] = 100;
							i--;
							j -= 2;
							if (i < 4)
								i = 4;
							if (j < 4)
								j = 4;
						}
					}
				}

				/* }}} */
			}
}

/* }}} */
/* {{{ susan_edges(in,r,sf,max_no,out) */

susan_edges(in, r, mid, bp, max_no, x_size, y_size)
	uchar *in, *bp, *mid;int *r, max_no, x_size, y_size; {
	float z;
	int do_symmetry, i, j, m, n, a, b, x, y, w;
	uchar c, *p, *cp;

	memset(r, 0, x_size * y_size * sizeof(int));

	//TODO
	for (i = 3; i < y_size - 3; i++)
		for (j = 3; j < x_size - 3; j++) {
			n = 100;
			p = in + (i - 3) * x_size + j - 1;
			cp = bp + in[i * x_size + j];

			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p);
			p += x_size - 3;

			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p);
			p += x_size - 5;

			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p);
			p += x_size - 6;

			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p);
			p += 2;
			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p);
			p += x_size - 6;

			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p);
			p += x_size - 5;

			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p);
			p += x_size - 3;

			n += *(cp - *p++);
			n += *(cp - *p++);
			n += *(cp - *p);

			if (n <= max_no)
				r[i * x_size + j] = max_no - n;
		}
	

	//TODO
	for (i = 4; i < y_size - 4; i++)
		for (j = 4; j < x_size - 4; j++) {
			if (r[i * x_size + j] > 0) {
				m = r[i * x_size + j];
				n = max_no - m;
				cp = bp + in[i * x_size + j];

				if (n > 600) {
					p = in + (i - 3) * x_size + j - 1;
					x = 0;
					y = 0;

					c = *(cp - *p++);
					x -= c;
					y -= 3 * c;
					c = *(cp - *p++);
					y -= 3 * c;
					c = *(cp - *p);
					x += c;
					y -= 3 * c;
					p += x_size - 3;

					c = *(cp - *p++);
					x -= 2 * c;
					y -= 2 * c;
					c = *(cp - *p++);
					x -= c;
					y -= 2 * c;
					c = *(cp - *p++);
					y -= 2 * c;
					c = *(cp - *p++);
					x += c;
					y -= 2 * c;
					c = *(cp - *p);
					x += 2 * c;
					y -= 2 * c;
					p += x_size - 5;

					c = *(cp - *p++);
					x -= 3 * c;
					y -= c;
					c = *(cp - *p++);
					x -= 2 * c;
					y -= c;
					c = *(cp - *p++);
					x -= c;
					y -= c;
					c = *(cp - *p++);
					y -= c;
					c = *(cp - *p++);
					x += c;
					y -= c;
					c = *(cp - *p++);
					x += 2 * c;
					y -= c;
					c = *(cp - *p);
					x += 3 * c;
					y -= c;
					p += x_size - 6;

					c = *(cp - *p++);
					x -= 3 * c;
					c = *(cp - *p++);
					x -= 2 * c;
					c = *(cp - *p);
					x -= c;
					p += 2;
					c = *(cp - *p++);
					x += c;
					c = *(cp - *p++);
					x += 2 * c;
					c = *(cp - *p);
					x += 3 * c;
					p += x_size - 6;

					c = *(cp - *p++);
					x -= 3 * c;
					y += c;
					c = *(cp - *p++);
					x -= 2 * c;
					y += c;
					c = *(cp - *p++);
					x -= c;
					y += c;
					c = *(cp - *p++);
					y += c;
					c = *(cp - *p++);
					x += c;
					y += c;
					c = *(cp - *p++);
					x += 2 * c;
					y += c;
					c = *(cp - *p);
					x += 3 * c;
					y += c;
					p += x_size - 5;

					c = *(cp - *p++);
					x -= 2 * c;
					y += 2 * c;
					c = *(cp - *p++);
					x -= c;
					y += 2 * c;
					c = *(cp - *p++);
					y += 2 * c;
					c = *(cp - *p++);
					x += c;
					y += 2 * c;
					c = *(cp - *p);
					x += 2 * c;
					y += 2 * c;
					p += x_size - 3;

					c = *(cp - *p++);
					x -= c;
					y += 3 * c;
					c = *(cp - *p++);
					y += 3 * c;
					c = *(cp - *p);
					x += c;
					y += 3 * c;

					z = sqrt((float) ((x * x) + (y * y)));
					if (z > (0.9 * (float) n)) /* 0.5 */
					{
						do_symmetry = 0;
						if (x == 0)
							z = 1000000.0;
						else
							z = ((float) y) / ((float) x);
						if (z < 0) {
							z = -z;
							w = -1;
						} else
							w = 1;
						if (z < 0.5) { /* vert_edge */
							a = 0;
							b = 1;
						} else {
							if (z > 2.0) { /* hor_edge */
								a = 1;
								b = 0;
							} else { /* diag_edge */
								if (w > 0) {
									a = 1;
									b = 1;
								} else {
									a = -1;
									b = 1;
								}
							}
						}
						if ((m > r[(i + a) * x_size + j + b])
								&& (m >= r[(i - a) * x_size + j - b])
								&& (m > r[(i + (2 * a)) * x_size + j + (2 * b)])
								&& (m >= r[(i - (2 * a)) * x_size + j - (2 * b)]))
							mid[i * x_size + j] = 1;
					} else
						do_symmetry = 1;
				} else
					do_symmetry = 1;

				if (do_symmetry == 1) {
					p = in + (i - 3) * x_size + j - 1;
					x = 0;
					y = 0;
					w = 0;

					/*   |      \
               y  -x-  w
					 |        \   */

					c = *(cp - *p++);
					x += c;
					y += 9 * c;
					w += 3 * c;
					c = *(cp - *p++);
					y += 9 * c;
					c = *(cp - *p);
					x += c;
					y += 9 * c;
					w -= 3 * c;
					p += x_size - 3;

					c = *(cp - *p++);
					x += 4 * c;
					y += 4 * c;
					w += 4 * c;
					c = *(cp - *p++);
					x += c;
					y += 4 * c;
					w += 2 * c;
					c = *(cp - *p++);
					y += 4 * c;
					c = *(cp - *p++);
					x += c;
					y += 4 * c;
					w -= 2 * c;
					c = *(cp - *p);
					x += 4 * c;
					y += 4 * c;
					w -= 4 * c;
					p += x_size - 5;

					c = *(cp - *p++);
					x += 9 * c;
					y += c;
					w += 3 * c;
					c = *(cp - *p++);
					x += 4 * c;
					y += c;
					w += 2 * c;
					c = *(cp - *p++);
					x += c;
					y += c;
					w += c;
					c = *(cp - *p++);
					y += c;
					c = *(cp - *p++);
					x += c;
					y += c;
					w -= c;
					c = *(cp - *p++);
					x += 4 * c;
					y += c;
					w -= 2 * c;
					c = *(cp - *p);
					x += 9 * c;
					y += c;
					w -= 3 * c;
					p += x_size - 6;

					c = *(cp - *p++);
					x += 9 * c;
					c = *(cp - *p++);
					x += 4 * c;
					c = *(cp - *p);
					x += c;
					p += 2;
					c = *(cp - *p++);
					x += c;
					c = *(cp - *p++);
					x += 4 * c;
					c = *(cp - *p);
					x += 9 * c;
					p += x_size - 6;

					c = *(cp - *p++);
					x += 9 * c;
					y += c;
					w -= 3 * c;
					c = *(cp - *p++);
					x += 4 * c;
					y += c;
					w -= 2 * c;
					c = *(cp - *p++);
					x += c;
					y += c;
					w -= c;
					c = *(cp - *p++);
					y += c;
					c = *(cp - *p++);
					x += c;
					y += c;
					w += c;
					c = *(cp - *p++);
					x += 4 * c;
					y += c;
					w += 2 * c;
					c = *(cp - *p);
					x += 9 * c;
					y += c;
					w += 3 * c;
					p += x_size - 5;

					c = *(cp - *p++);
					x += 4 * c;
					y += 4 * c;
					w -= 4 * c;
					c = *(cp - *p++);
					x += c;
					y += 4 * c;
					w -= 2 * c;
					c = *(cp - *p++);
					y += 4 * c;
					c = *(cp - *p++);
					x += c;
					y += 4 * c;
					w += 2 * c;
					c = *(cp - *p);
					x += 4 * c;
					y += 4 * c;
					w += 4 * c;
					p += x_size - 3;

					c = *(cp - *p++);
					x += c;
					y += 9 * c;
					w -= 3 * c;
					c = *(cp - *p++);
					y += 9 * c;
					c = *(cp - *p);
					x += c;
					y += 9 * c;
					w += 3 * c;

					if (y == 0)
						z = 1000000.0;
					else
						z = ((float) x) / ((float) y);
					if (z < 0.5) { /* vertical */
						a = 0;
						b = 1;
					} else {
						if (z > 2.0) { /* horizontal */
							a = 1;
							b = 0;
						} else { /* diagonal */
							if (w > 0) {
								a = -1;
								b = 1;
							} else {
								a = 1;
								b = 1;
							}
						}
					}
					if ((m > r[(i + a) * x_size + j + b])
							&& (m >= r[(i - a) * x_size + j - b])
							&& (m > r[(i + (2 * a)) * x_size + j + (2 * b)])
							&& (m >= r[(i - (2 * a)) * x_size + j - (2 * b)]))
						mid[i * x_size + j] = 2;
				}
			}
		}
}

/* }}} */

/* {{{ main(argc, argv) */

main(argc, argv)
	int argc;char *argv[]; {

	uchar *in, *bp, *mid;
	uchar **image;
	int *r, bt = 20, drawing_mode = 0, max_no_edges = 2650, mode = 0, n_cols, n_rows;

	int    proc_id, n_procs, n_workers, rows_per_worker, offset, source, destination, start, end, i, j, h;
 
  	MPI_Status status;       /* MPI structure containing return codes
                              for message passing operations         */
	
	if (argc < 3)
		printf("Usage: susan <in.pgm> <out.pgm>\n\n");

	MPI_Init(&argc, &argv); 
  
  	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
 
  	n_workers = n_procs - 1;

	if (proc_id == MASTER) {

		get_image(argv[1], &in, &n_cols, &n_rows);

		rows_per_worker = (n_rows % n_workers == 0) ? (n_rows / n_workers) : (n_rows / n_workers + 1);

		offset = 0;

		for (i = 1; i <= n_workers; i++) {
      		destination = i;

			if(n_rows % n_workers != 0 && i == n_workers)
				 rows_per_worker = n_rows % n_workers;
      
      		MPI_Send(&n_cols,           1, 							MPI_INT, 			destination, BEGIN, MPI_COMM_WORLD);
      		MPI_Send(&rows_per_worker,  1, 							MPI_INT, 			destination, BEGIN, MPI_COMM_WORLD); 
      		MPI_Send(in + offset, 		rows_per_worker * n_cols, 	MPI_UNSIGNED_CHAR, 	destination, BEGIN, MPI_COMM_WORLD);
 
      		offset += rows_per_worker * n_cols;
    	}

		rows_per_worker = (n_rows % n_workers == 0) ? (n_rows / n_workers) : (n_rows / n_workers + 1);

		for(i = 1; i <= n_workers; i++) {         
        	source = i;

			offset = (i - 1) * rows_per_worker * n_cols;

			if(n_rows % n_workers != 0 && i == n_workers)
				 rows_per_worker = n_rows % n_workers;
 
        	MPI_Recv(in + offset, rows_per_worker * n_cols, MPI_UNSIGNED_CHAR, source, DONE, MPI_COMM_WORLD, &status);
      	}

		put_image(argv[2], in, n_cols, n_rows);

	}

	else {

		MPI_Recv(&n_cols,	1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
    	MPI_Recv(&n_rows,  	1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);

		in = (uchar *) malloc(n_cols * n_rows);

		MPI_Recv(in, n_rows * n_cols, MPI_UNSIGNED_CHAR, MASTER, BEGIN, MPI_COMM_WORLD, &status);

		r = (int *) malloc(n_cols * n_rows * sizeof(int));
		setup_brightness_lut(&bp, bt, 6);

		mid = (uchar *) malloc(n_cols * n_rows);
		memset(mid, 100, n_cols * n_rows); /* note not set to zero */

		susan_edges(in, r, mid, bp, max_no_edges, n_cols, n_rows);

		susan_thin(r, mid, n_cols, n_rows);
		edge_draw(in, mid, n_cols, n_rows, drawing_mode);
		
		MPI_Send(in, n_rows * n_cols, MPI_UNSIGNED_CHAR, MASTER, DONE, MPI_COMM_WORLD);
	}

	MPI_Finalize();
}
