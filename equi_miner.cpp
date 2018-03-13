// Wagner's algorithm for Generalized Birthday Paradox, a memory-hard proof-of-work
// Copyright (c) 2016 John Tromp

#include "equi_miner.h"
#include <unistd.h>
#include "ctype.h"

int hextobyte(const char * x) {
  u32 b = 0;
  for (int i = 0; i < 2; i++) {
    uchar c = tolower(x[i]);
    assert(isxdigit(c));
    b = (b << 4) | (c - (c >= '0' && c <= '9' ? '0' : ('a' - 10)));
  }
  return b;
}

/*
** Write ZCASH_SOL_LEN bytes representing the encoded solution as per the
** Zcash protocol specs (512 x 21-bit inputs).
**
** out		ZCASH_SOL_LEN-byte buffer where the solution will be stored
** inputs	array of 32-bit inputs
** n		number of elements in array
*/
#define PREFIX (DIGITBITS)
void store_encoded_sol(uint8_t *out, uint32_t *inputs, uint32_t n)
{
    uint32_t byte_pos = 0;
    int32_t bits_left = PREFIX + 1;
    uint8_t x = 0;
    uint8_t x_bits_used = 0;
    while (byte_pos < n)
      {
        if (bits_left >= 8 - x_bits_used)
          {
            x |= inputs[byte_pos] >> (bits_left - 8 + x_bits_used);
            bits_left -= 8 - x_bits_used;
            x_bits_used = 8;
          }
        else if (bits_left > 0)
          {
            uint32_t mask = ~(-1 << (8 - x_bits_used));
            mask = ((~mask) >> bits_left) & mask;
            x |= (inputs[byte_pos] << (8 - x_bits_used - bits_left)) & mask;
            x_bits_used += bits_left;
            bits_left = 0;
          }
        else if (bits_left <= 0)
          {
            assert(!bits_left);
            byte_pos++;
            bits_left = PREFIX + 1;
          }
        if (x_bits_used == 8)
          {
	    *out++ = x;
            x = x_bits_used = 0;
          }
      }
}

#define ZCASH_SOL_LEN                   ((1 << WK) * (PREFIX + 1) / 8)
void print_encoded_sol(uint32_t *inputs, uint32_t n)
{
    uint8_t	sol[ZCASH_SOL_LEN];
    uint32_t	i;
    store_encoded_sol(sol, inputs, n);
    for (i = 0; i < sizeof (sol); i++)
	printf("%02x", sol[i]);
    printf("\n");
    fflush(stdout);
}

int main(int argc, char **argv) {
  int nthreads = 1;
  int nonce = 0;
  int range = 1;
  bool showsol = false;
  const char *header = "";
  const char *hex = "";
  int noncepos = 32;
  int c;
  while ((c = getopt (argc, argv, "p:h:n:r:t:x:s")) != -1) {
    switch (c) {
      case 'h':
        header = optarg;
        printf("Header specified: %s\n", header);
        break;
      case 'n':
        nonce = atoi(optarg);
        printf("Nonce specified: %d\n", nonce);
        break;
      case 'r':
        range = atoi(optarg);
        printf("Range specified: %d\n", range);
        break;
      case 's':
        showsol = true;
        break;
      case 't':
        nthreads = atoi(optarg);
        break;
      case 'x':
        hex = optarg;
        printf("Hex specified: %s\n", hex);
        break;
      case 'p':
        noncepos = atoi(optarg);
        break;
    }
  }
  printf("NoncePos: %d\n", noncepos);
#ifndef XWITHASH
  if (sizeof(tree) > 4)
    printf("WARNING: please compile with -DXWITHASH to shrink tree!\n");
#endif
#ifdef ATOMIC
  if (nthreads==1)
    printf("WARNING: use of atomics hurts single threaded performance!\n");
#else
  assert(nthreads==1);
#endif
  printf("Looking for wagner-tree on (\"%s\",%d", hex ? "0x..." : header, nonce);
  if (range > 1)
    printf("-%d", nonce+range-1);
  printf(") with %d %d-bit digits and %d threads\n", NDIGITS, DIGITBITS, nthreads);
  thread_ctx *threads = (thread_ctx *)calloc(nthreads, sizeof(thread_ctx));
  assert(threads);
  equi eq(nthreads);
  printf("Using 2^%d buckets, %dMB of memory, and %d-way blake2b\n", BUCKBITS, 1 + eq.hta.alloced / 0x100000, NBLAKES);
#ifdef ASM_BLAKE
  printf("Using xenoncat's assembly blake code\n");
#endif
  u32 sumnsols = 0;
  char headernonce[HEADERNONCELEN];
  u32 hdrlen = strlen(header);
  if (*hex) {
    assert(strlen(hex) == 2 * HEADERNONCELEN);
    for (int i = 0; i < HEADERNONCELEN; i++)
      headernonce[i] = hextobyte(&hex[2*i]);
  } else {
    memcpy(headernonce, header, hdrlen);
    memset(headernonce+hdrlen, 0, sizeof(headernonce)-hdrlen);
  }

  for (int r = 0; r < range; r++) {
    ((u32 *)headernonce)[noncepos] = htole32(nonce+r);
    eq.setheadernonce(headernonce, sizeof(headernonce));

    for (int t = 0; t < nthreads; t++) {
      threads[t].id = t;
      threads[t].eq = &eq;
      int err = pthread_create(&threads[t].thread, NULL, worker, (void *)&threads[t]);
      assert(err == 0);
    }
    for (int t = 0; t < nthreads; t++) {
      int err = pthread_join(threads[t].thread, NULL);
      assert(err == 0);
    }
    u32 nsols, maxsols = min(MAXSOLS, eq.nsols);
    for (nsols = 0; nsols < maxsols; nsols++) {
      if (showsol) {
        printf("\nSolution [%d]\n", PROOFSIZE);
        for (u32 i = 0; i < PROOFSIZE; i++)
          printf(" %jx", (uintmax_t)eq.sols[nsols][i]);

        printf("\nEncoded\n");
        print_encoded_sol(eq.sols[nsols], 1 << WK);
      }
    }
    if (maxsols > 0)
    {
      printf("\nNonce ");
      for (int i = 0; i < 32; ++i)
      {
        printf("%02x", (unsigned char)headernonce[108+i]);
      }
      printf(" %d solutions\n", nsols);
    }
    sumnsols += nsols;
  }
  free(threads);
  printf("%d total solutions\n", sumnsols);
  return 0;
}
