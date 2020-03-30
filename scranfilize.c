// Copyright (C) 2018-2020, Armin Biere, Johannes Kepler University Linz, Austria

const char * usage =
"usage: scranfilize [ <option> ... ] [ <original-cnf> [ <scrambled-cnf> ] ]\n"
"\n"
"where '<option>' is one of the following\n"
"\n"
"   -h         print this command line option summary and exit\n"
"   --version  print version and exit\n"
"\n"
"   -p         completely permute variables\n"
"   -P         completely permute clauses\n"
"\n"
"   -r         reverse order of all variables\n"
"   -R         reverse order of all clauses\n"
"\n"
"   -s <seed>  random number generator seed\n"
"              (default is to hash time and process id)\n"
"\n"
"   -f <prob>  probability of flipping a literal (default '.01')\n"
"   -v <win>   relative variable move window (default '.01')\n"
"   -c <win>   relative clause move window (default '.01')\n"
"\n"
"   -a         use absolute move windows (window defaults become '1')\n"
"\n"
"   --force    force to overwrite existing file\n"
"\n"
"by default the original CNF is read from '<stdin>' unless '<original-cnf>'\n"
"is given.  The scrambled CNF is written to '<stdout>' or '<scrambled-cnf>'.\n"
;

/*------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/times.h>
#include <sys/types.h>
#include <unistd.h>

/*------------------------------------------------------------------------*/

#include "config.h"

/*------------------------------------------------------------------------*/

// Options.

static long seed = -1;
static bool permute_variables = false;
static bool permute_clauses = false;
static bool reverse_variables = false;
static bool reverse_clauses = false;
static double literal_flip_probability = -1;
static double variable_move_window = -1;
static double clause_move_window = -1;
static bool absolute_windows = false;
static bool force = false;

/*------------------------------------------------------------------------*/

// CNF.

static int max_var;
static int num_clauses;
static int ** clauses;

/*------------------------------------------------------------------------*/

// Scrambling maps.

static bool * flipped;
static int * clause_map;
static int * variable_map;

/*------------------------------------------------------------------------*/

static void die (const char * msg, ...) {
  fflush (stdout);
  fputs ("scranfilize: error: ", stderr);
  va_list ap;
  va_start (ap, msg);
  vfprintf (stderr, msg, ap);
  va_end (ap);
  fputc ('\n', stderr);
  exit (1);
}

static void msg (const char * msg, ...) {
  fflush (stdout);
  fputs ("[scranfilize] ", stderr);
  va_list ap;
  va_start (ap, msg);
  vfprintf (stderr, msg, ap);
  va_end (ap);
  fputc ('\n', stderr);
  fflush (stderr);
}

/*------------------------------------------------------------------------*/

bool exists_file (const char * path) {
  struct stat buf;
  return !stat (path, &buf);
}

static bool is_suffix (const char * path, const char * suffix) {
  size_t k = strlen (path), l = strlen (suffix);
  return k > l && !strcmp (path + k - l, suffix);
}

static FILE *
open_pipe (const char * path, const char * fmt, int * close_file) {
  char * cmd = malloc (strlen (fmt) + strlen (path));
  if (!cmd) die ("out-of-memory allocating command string");
  sprintf (cmd, fmt, path);
  FILE * file = popen (cmd, "r");
  *close_file = 2;
  free (cmd);
  return file;
}

static int next_char (FILE * file, int * lineno) {
  int res = getc (file);
  if (res == '\n') *lineno += 1;
  return res;
}

static bool space (int ch) {
  return ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n';
}

static void
parse_error (const char * path, int lineno, const char * msg, ...) {
  fflush (stdout);
  fprintf (stderr,
    "scranfilize: parse error: %s:%d: ",
    path, lineno);
  va_list ap;
  va_start (ap, msg);
  vfprintf (stderr, msg, ap);
  va_end (ap);
  fputc ('\n', stderr);
  exit (1);
}

static void parse (const char * path) {

#define suffix(STR) is_suffix (path, STR)
#define pipe(CMD) file = open_pipe (path, CMD, &close_file)
#define next() next_char (file, &lineno)
#define perr(...) parse_error (path, lineno, __VA_ARGS__)

  if (path && !exists_file (path)) die ("file '%s' does not exist", path);

  FILE * file;
  int close_file;

  if (!path) {
    path = "<stdin>";
    file = stdin;
    close_file = 0;
  } else if (suffix (".xz") || suffix (".lzma")) pipe ("xz -c -d %s");
  else if (suffix (".bz2")) pipe ("bzip2 -c -d %s");
  else if (suffix (".gz")) pipe ("gzip -c -d %s");
  else if (suffix (".7z")) pipe ("7z x -so %s 2>/dev/null");
  else {
    file = fopen (path, "r");
    close_file = 1;
  }
  if (!file) die ("can not read original CNF '%s'", path);
  msg ("reading original CNF from '%s'", path);

  int lineno = 1;

  int ch;

  for (;;) {
    ch = next ();
    if (ch == EOF) perr ("unexpected end-of-file before header");
    if (ch == 'p') break;
    if (ch == 'c') {
      while ((ch = next ()) != '\n')
	if (ch == EOF)
	  perr ("unexpected end-of-file in header comment");
      continue;
    }
    if (ch == EOF) perr ("unexpected end-of-file");
    else if (isprint (ch)) perr ("unexpected character '%c'", ch);
    else perr ("unexpected character (code '%d')", ch);
  }

  assert (ch == 'p');
  if (next () != ' ' ||
      next () != 'c' ||
      next () != 'n' ||
      next () != 'f' ||
      next () != ' ')
    perr ("invalid DIMACS header");

  ch = next ();
  if (!isdigit (ch)) perr ("expected digit after 'p cnf '");
  max_var = ch - '0';
  while (isdigit (ch = next ())) {
    if (INT_MAX/10 < max_var) perr ("variable number way too large");
    max_var *= 10;
    const int digit = ch - '0';
    if (INT_MAX - digit < max_var) perr ("variable number too large");
    max_var += digit;
  }

  if (ch != ' ') perr ("expected space after variable number");

  ch = next ();
  if (!isdigit (ch)) perr ("expected digit after 'p cnf %d'", max_var);

  int specified_clauses = ch - '0';
  while (isdigit (ch = next ())) {
    if (INT_MAX/10 < specified_clauses)
      perr ("clause number way too large");
    specified_clauses *= 10;
    const int digit = ch - '0';
    if (INT_MAX - digit < specified_clauses)
      perr ("clause number too large");
    specified_clauses += digit;
  }

  msg ("found 'p cnf %d %d' header", max_var, specified_clauses);

  while (ch != '\n') {
    if (!space (ch)) perr ("expected white space before new line");
    ch = next ();
  }

  clauses = malloc (specified_clauses * sizeof *clauses);
  if (!clauses) die ("out-of-memory allocating clauses");

  int num_literals = 0, size_literals = 0, * literals = 0;

  ch = next ();
  for (;;) {
    if (space (ch)) ch = next ();
    else if (ch == EOF) {
      if (num_literals) perr ("terminating zero missing");
      if (num_clauses < specified_clauses)
	perr ("%d clause%s missing",
	  num_clauses,
	  num_clauses + 1 == specified_clauses ? "" : "s");
      break;
    } else if (ch == 'c') {
      while ((ch = next ()) != '\n' && ch != EOF)
	;
    } else {
      int sign;
      if (ch == '-') {
	ch = next ();
	if (!isdigit (ch)) perr ("expected digit after '-'");
	if (ch == '0') perr ("expected non-zer digit after '-'");
	sign = -1;
      } else {
	if (!isdigit (ch)) perr ("expected digit or '-'");
	sign = 1;
      }
      assert (isdigit (ch));
      int idx = ch - '0';
      while (isdigit (ch = next ())) {
	if (INT_MAX/10 < idx)
	  perr ("variable way too large");
	idx *= 10;
	const int digit = ch - '0';
	if (INT_MAX - digit < idx)
	  perr ("variable too large");
	idx += digit;
      }
      if (idx > max_var) perr ("maximum variable index exceeded");
      if (!space (ch) && ch != 'c' && ch != EOF) {
	if (isprint (ch))
	  perr ("unexpected character '%c' after literal", ch);
	else
	  perr ("unexpected character after literal (code '%d')", ch);
      }
      if (num_clauses == specified_clauses) perr ("too many clauses");
      int lit = sign * idx;
      if (lit) {
	if (num_literals == size_literals) {
	  if (size_literals) size_literals *= 2; else size_literals = 1;
	  literals = realloc (literals, size_literals * sizeof *literals);
	  if (!literals)
	    die ("out-of-memory reallocating literal stack");
	}
	literals[num_literals++] = lit;
      } else {
	int * clause = malloc ((num_literals + 1) * sizeof *clause);
	if (!clause) die ("out-of-memory allocating clause");
	for (int i = 0; i < num_literals; i++) clause[i] = literals[i];
	clause[num_literals] = 0;
	num_literals = 0;
	assert (num_clauses < specified_clauses);
	clauses[num_clauses++] = clause;
      }
    }
  }

  if (close_file == 1) fclose (file);
  if (close_file == 2) pclose (file);
  if (literals) free (literals);
}

/*------------------------------------------------------------------------*/

typedef struct Map { int src; double dst; } Map;

static int cmp_rank (const void * p, const void * q) {
  Map * r = (Map *) p, * s = (Map *) q;
  if (r->dst < s->dst) return -1;
  if (r->dst > s->dst) return 1;
  if (r->src < s->src) return -1;
  if (r->src > s->src) return 1;
  return 0;
}

static int * rank (int n, bool permute, double width) {

  srand48 (seed);

  typedef struct Map { int src; double dst; } Map;

  Map * ranks = malloc (n * sizeof *ranks);
  if (!ranks) die ("out-of-memory allocating %d ranks", n);

  for (int i = 0; i < n; i++) {
    Map * m = ranks + i;
    m->src = i;
    if (permute) m->dst = drand48 () * n;
    else {
      double tmp = drand48 () * width;
      if (!absolute_windows) tmp *= n;
      m->dst = i + tmp;
    }
  }

  qsort (ranks, n, sizeof *ranks, cmp_rank);

#if 0
do {
  static int print = 0;
  if (print++) break;
  for (int i = 0; i < n; i++) {
    Map * m = ranks + i;
    printf ("%d %f\n", m->src, m->dst);
  }
} while (0);
#endif

  int * res = malloc (n * sizeof *res);
  if (!res) die ("out-of-memory allocating %d map", n);

  for (int i = 0; i < n; i++)
    res[i] = ranks[i].src;

  free (ranks);

#if 0
do {
  static int print = 0;
  if (print++) break;
  for (int i = 0; i < n; i++)
    printf ("%d %d\n", i, res[i]);
} while (0);
#endif


#if 0
  {
#   define DATA "/tmp/scranfilize-data"
#   define CMD "/tmp/scranfilize-CMD"
    FILE * file = fopen (DATA, "w");
    for (int i = 0; i < n; i++) fprintf (file, "%d %d\n", i, res[i]);
    fclose (file);
    file = fopen (CMD, "w");
    fprintf (file, "plot \"" DATA "\"\npause mouse\nquit\n");
    fclose (file);
    system ("gnuplot " CMD);
  }
#endif

  return res;
}

/*------------------------------------------------------------------------*/

static bool * flip () {
  srand48 (seed);

  bool * res = malloc (max_var * sizeof *res);
       if (literal_flip_probability <= 0.0) memset (res, 0, max_var);
  else if (literal_flip_probability >= 1.0) memset (res, 1, max_var);
  else {
    srand48 (seed);
    for (int i = 0; i < max_var; i++)
      res[i] = (drand48 () <= literal_flip_probability);
  }

#if 0
  for (int i = 0; i < max_var; i++)
    printf ("%d %d\n", i, (int) res[i]);
#endif

  return res;
}

/*------------------------------------------------------------------------*/

static void scramble () {
  variable_map = rank (max_var, permute_variables, variable_move_window);
  clause_map = rank (num_clauses, permute_clauses, clause_move_window);
  flipped = flip ();
}

/*------------------------------------------------------------------------*/

static void banner (FILE * file, 
                    void (*print)(FILE * file, const char *, ...)) {
  print (file, "Scranfilize CNF Scrambler");
  print (file, "Version %s %s", VERSION, GITID);
  print (file, "random seed '%ld'", seed);
  if (reverse_variables) print (file, "reverse all clauses ('-r')");
  if (reverse_clauses) print (file, "reverse all variables ('-R')");
  print (file, "literal flip probability %g ('-f %g')",
    literal_flip_probability, literal_flip_probability);
  if (permute_variables)
    print (file, "randomly permuting variables");
  else
    print (file, "%s variable move window %g ('-v %g')",
      absolute_windows ? "absolute" : "relative",
      variable_move_window, variable_move_window);
  if (permute_clauses)
    print (file, "randomly permuting clauses");
  else
    print (file, "%s clause move window %g ('-c %g')",
      absolute_windows ? "absolute" : "relative",
      clause_move_window, clause_move_window);
}

/*------------------------------------------------------------------------*/

static bool exists (const char * path) {
  struct stat buf;
  return !stat (path, &buf);
}

static void print_message (FILE *file, const char * msg, ...) {
  fputs ("c ", file);
  va_list ap;
  va_start (ap, msg);
  vfprintf (file, msg, ap);
  va_end (ap);
  fputc ('\n', file);
}

static void print (const char * path) {

  if (path && exists (path)) {
    if (force) msg ("forced to overwrite existing '%s'", path);
    else die ("path '%s' exist (use '--force')", path);
  }

  FILE * file;
  int close_file;

  if (path) {
    if (!(file = fopen (path, "w")))
      die ("can not write scrambled CNF '%s'", path);
    close_file = 1;
  } else {
    path = "<stdout>";
    file = stdout;
    close_file = 0;
  }

  msg ("writing scrambled CNF to '%s'", path ? path : "<stdout>");

  banner (file, print_message);

  fprintf (file, "p cnf %d %d\n", max_var, num_clauses);
  for (int i = 0; i < num_clauses; i++) {
    int j = clause_map[i];
    if (reverse_clauses) j = num_clauses-1 - j;
    assert (0 <= j), assert (j < num_clauses);
    for (const int * p = clauses[j]; *p; p++) {
      const int src = *p;
      int idx = abs (src);
      if (reverse_variables) idx = max_var + 1 - idx;
      assert (1 <= idx), assert (idx <= max_var);
      int dst = variable_map[idx-1] + 1;
      assert (1 <= dst), assert (dst <= max_var);
      if (src < 0) dst = -dst;
      if (flipped[idx-1]) dst = -dst;
      fprintf (file, "%d ", dst);
    }
    fprintf (file, "0\n");
  }

  if (close_file == 1) fclose (file);
  if (close_file == 2) pclose (file);
}

/*------------------------------------------------------------------------*/
// Files.

static const char * original;
static const char * scrambled;

/*------------------------------------------------------------------------*/

static bool valid (double f) {
  if (f < 0 && (f < -1e150 || f > -1e-150))
    return false;
  if (f > 0 && (f > +1e150 || f < 1e-150))
    return false;
  return true;
}

static void init (int argc, char ** argv) {

  for (int i = 1; i < argc; i++) {
    if (!strcmp (argv[i], "-h")) fputs (usage, stdout), exit (0);
    else if (!strcmp (argv[i], "--version"))
      printf ("%s\n", VERSION), exit (0);
    else if (!strcmp (argv[i], "-p")) permute_variables = true;
    else if (!strcmp (argv[i], "-P")) permute_clauses = true;
    else if (!strcmp (argv[i], "-r")) reverse_variables = true;
    else if (!strcmp (argv[i], "-R")) reverse_clauses = true;
    else if (!strcmp (argv[i], "-s")) {
      if (++i == argc) die ("argument to '-s' missing");
      seed = atol (argv[i]);
      if (seed < 0) die ("invalid negative argument to '-s'");
    } else if (!strcmp (argv[i], "-f")) {
      if (++i == argc) die ("argument to '-f' missing");
      double tmp = atof (argv[i]);
      if (!valid (tmp) || tmp > 1.0)
	die ("invalid argument in '-f %s'", argv[i]);
      if (literal_flip_probability >= 0)
	die ("multiple '-f' options");
      literal_flip_probability = tmp;
    } else if (!strcmp (argv[i], "-v")) {
      if (++i == argc) die ("argument to '-v' missing");
      double tmp = atof (argv[i]);
      if (!valid (tmp))
	die ("invalid argument in '-v %s'", argv[i]);
      if (variable_move_window >= 0)
	die ("multiple '-v' options");
      variable_move_window = tmp;
    } else if (!strcmp (argv[i], "-c")) {
      if (++i == argc) die ("argument to '-c' missing");
      double tmp = atof (argv[i]);
      if (!valid (tmp))
	die ("invalid argument in '-c %s'", argv[i]);
      if (clause_move_window >= 0)
	die ("multiple '-c' options");
      clause_move_window = tmp;
    } else if (!strcmp (argv[i], "-a")) absolute_windows = true;
    else if (!strcmp (argv[i], "--force")) force = true;
    else if (argv[i][0] == '-')
      die ("invalid option '%s' (try '-h')", argv[i]);
    else if (scrambled)
      die ("too many arguments '%s', '%s' and '%s' (try '-h')",
        original, scrambled, argv[i]);
    else if (original) scrambled = argv[i];
    else original = argv[i];
  }

  if (permute_variables) {
    if (reverse_variables) die ("can not combine '-p' and '-r'");
    if (variable_move_window >=0) die ("can not combine '-p' and '-v'");
    if (absolute_windows) die ("can not combine '-p' and '-a'");
  }

  if (permute_clauses) {
    if (reverse_clauses) die ("can not combine '-P' and '-R'");
    if (clause_move_window >=0) die ("can not combine '-P' and '-c'");
    if (absolute_windows) die ("can not combine '-P' and '-a'");
  }

  if (seed < 0) {
    uint64_t t = 8526563 * (unsigned long) times (0);
    uint64_t p = 3944621 * (unsigned long) getpid ();
    uint64_t tmp = t + p;
    seed = tmp & 0xffffffff;
    seed ^= tmp >> 32;
  }

  if (literal_flip_probability < 0) literal_flip_probability = 0.01;

  const double default_window = absolute_windows ? 1.0 : 0.01;
  if (variable_move_window < 0) variable_move_window = default_window;
  if (clause_move_window < 0) clause_move_window = default_window;

  banner (stdout, print_message);
}

/*------------------------------------------------------------------------*/

void reset () {
  free (flipped);
  free (clause_map);
  free (variable_map);
  for (int i = 0; i < num_clauses; i++) free (clauses[i]);
  free (clauses);
}

/*------------------------------------------------------------------------*/

int main (int argc, char ** argv) {
  init (argc, argv);
  parse (original);
  scramble ();
  print (scrambled);
  reset ();
  return 0;
}
