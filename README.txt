CHALLENGE IMPLEMENTATION

The challenge is composed of 5 files.

1 - The file AFF_MULT_challenge.bin is composed of the polynomials representing the affine multiple of the challenge. It is represented as a n by n matrix of degree 3 polynomials. They are stored by lines. The polynomials are stored in deglex order. The polynomials are represented by their monomials, and grouped by 8 to be encoded tightly into 8 bits. To see more about how they are used for signature, see MAIN_challenge.iypnb.

2 - The file PK_challenge.bin contains the public key of the challenge. Same representation as for AFF_MULT_challenge.bin exccept there are only n-a polynomials.

3 - The file TEST_challenge.bin contains the test values to verify that the affine multiple is functionnal.

4 - The file MAIN_challenge.iypnb contains a sample code to manipulate the affine multiple and sign with it.

5 - The file lib.sage contains some functions to manipulate the the previous files and sign.
