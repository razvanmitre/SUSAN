/*#define FOPENB*//* uncomment if using djgpp gnu C for DOS or certain Win95 compilers */
#define SEVEN_SUPP           /* size for non-max corner suppression; SEVEN_SUPP or FIVE_SUPP */
#define MAX_CORNERS   15000  /* max corners per frame */


#include <omp.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/file.h>    /* may want to remove this line */
#include <malloc.h>      /* may want to remove this line */
#define  exit_error(IFB,IFC) { fprintf(stderr,IFB,IFC); exit(0); }
#define  FTOI(a) ( (a) < 0 ? ((int)(a-0.5)) : ((int)(a+0.5)) )
typedef unsigned char uchar;



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



void setup_brightness_lut(bp, thresh, form)
	uchar **bp;int thresh, form; {
	int k;
	float temp;

	*bp = (unsigned char *) malloc(516);
	*bp = *bp + 258;

	for (k = -256; k < 257; k++) {
		temp = ((float) k) / ((float) thresh);
		temp = temp * temp;
		if (form == 6)
			temp = temp * temp * temp;
		temp = 100.0 * exp(-temp);
		*(*bp + k) = (uchar) temp;
	}
}



edge_draw(in, mid, x_size, y_size, drawing_mode)
	uchar *in, *mid;int x_size, y_size, drawing_mode; {
	int i;
	uchar *inp, *midp;

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



susan_thin(r, mid, x_size, y_size)
	uchar *mid;int *r, x_size, y_size; {
	int l[9], centre, nlinks, npieces, b01, b12, b21, b10, p1, p2, p3, p4, b00,
			b02, b20, b22, m, n, a, b, x, y, i, j;
	uchar *mp;

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
					} else if (mid[(i - 1) * x_size + j] < 8) {
						l[1] = 0;
						l[0] = 0;
						l[2] = 0;
						l[3] *= 2;
						l[5] *= 2;
						l[6] *= 3;
						l[8] *= 3;
						l[7] *= 4;
					} else if (mid[(i - 1) * x_size + j + 1] < 8) {
						l[2] = 0;
						l[1] = 0;
						l[5] = 0;
						l[0] *= 2;
						l[8] *= 2;
						l[3] *= 3;
						l[7] *= 3;
						l[6] *= 4;
					} else if (mid[(i) * x_size + j - 1] < 8) {
						l[3] = 0;
						l[0] = 0;
						l[6] = 0;
						l[1] *= 2;
						l[7] *= 2;
						l[2] *= 3;
						l[8] *= 3;
						l[5] *= 4;
					} else if (mid[(i) * x_size + j + 1] < 8) {
						l[5] = 0;
						l[2] = 0;
						l[8] = 0;
						l[1] *= 2;
						l[7] *= 2;
						l[0] *= 3;
						l[6] *= 3;
						l[3] *= 4;
					} else if (mid[(i + 1) * x_size + j - 1] < 8) {
						l[6] = 0;
						l[3] = 0;
						l[7] = 0;
						l[0] *= 2;
						l[8] *= 2;
						l[1] *= 3;
						l[5] *= 3;
						l[2] *= 4;
					} else if (mid[(i + 1) * x_size + j] < 8) {
						l[7] = 0;
						l[6] = 0;
						l[8] = 0;
						l[3] *= 2;
						l[5] *= 2;
						l[0] *= 3;
						l[2] *= 3;
						l[1] *= 4;
					} else if (mid[(i + 1) * x_size + j + 1] < 8) {
						l[8] = 0;
						l[5] = 0;
						l[7] = 0;
						l[6] *= 2;
						l[2] *= 2;
						l[1] *= 3;
						l[3] *= 3;
						l[0] *= 4;
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



susan_edges(in, r, mid, bp, max_no, x_size, y_size)
	uchar *in, *bp, *mid;int *r, max_no, x_size, y_size; {
	float z;
	int do_symmetry, i, j, m, n, a, b, x, y, w;
	uchar c, *p, *cp;

	memset(r, 0, x_size * y_size * sizeof(int));

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



main(int argc, char *argv[])
{
	char filename[80];
	uchar *image, *bp;
	int max_no_edges = 2650, 
		x_size,	y_size;

	int thread_id, nthreads,
		x_chunk_size, y_chunk_size,
		*r_chunk,
		offset_up = 0, offset_down = 0;
	uchar *image_chunk, *mid_chunk;


	// nthreads = omp_get_num_threads();
	nthreads = 4;

	get_image(argv[1], &image, &x_size, &y_size);
	setup_brightness_lut(&bp, 20, 6);

	printf("%d %d\n", nthreads, x_chunk_size);


	#pragma omp parallel private(thread_id, r_chunk, image_chunk, mid_chunk, offset_up, offset_down, x_chunk_size, y_chunk_size) shared(nthreads, image, x_size, y_size, max_no_edges, bp)
	{
		thread_id = omp_get_thread_num();
		y_chunk_size = (y_size % nthreads == 0) ? (y_size / nthreads) : (y_size / nthreads + 1);
		x_chunk_size = x_size;

		if(thread_id == 0)
		{
			offset_up = 0;
			offset_down = 8;
		}
		else if(thread_id == nthreads - 1)
		{
			offset_up = 6;
			offset_down = 0;
		}
		else
		{
			offset_up = 6;
			offset_down = 8;	
		}
		printf("thread: %d\n", thread_id);



		image_chunk = (uchar *) malloc(x_chunk_size * (y_chunk_size + offset_up + offset_down));
		memcpy(image_chunk, image + thread_id * y_chunk_size * x_size - offset_up * x_chunk_size, y_chunk_size * x_size + offset_down * x_chunk_size);


		y_chunk_size = y_chunk_size + offset_up + offset_down;

		r_chunk = (int *) malloc(x_chunk_size * y_chunk_size * sizeof(int));
		mid_chunk = (uchar *) malloc(x_chunk_size * y_chunk_size);
		memset(mid_chunk, 100, x_chunk_size * y_chunk_size);

		susan_edges(image_chunk, r_chunk, mid_chunk, bp, max_no_edges, x_chunk_size, y_chunk_size);
		susan_thin(r_chunk, mid_chunk, x_chunk_size, y_chunk_size);
		edge_draw(image_chunk, mid_chunk, x_chunk_size, y_chunk_size, 0);

		y_chunk_size = y_chunk_size - offset_up - offset_down;


		memcpy(image + thread_id * y_chunk_size * x_size, image_chunk + offset_up * x_chunk_size, y_chunk_size * x_size);
	}


	put_image(argv[2], image, x_size, y_size);
}


