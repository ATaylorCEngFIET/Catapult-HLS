#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "sbox.h"
#include "aes_enc.h"

void shift_row_enc(uint8_t state[4][4], uint8_t result[4][4]);
void subbytes(uint8_t state[4][4], uint8_t result[4][4]);
void addroundkey(uint8_t state[4][4],uint8_t iteration, uint8_t result[4][4],uint8_t ekey[240]);
void mixcolumn(uint8_t state[4][4], uint8_t result[4][4]);


void aes_enc(uint8_t state[16],uint8_t cipher[16],uint8_t ekey[240])

{

	uint8_t iteration = 0;
	uint8_t x,y;
	uint8_t sub[4][4];
	uint8_t shift[4][4];
	uint8_t mix[4][4];
	uint8_t round[4][4];
	uint8_t state_grid[4][4];
	uint8_t result[4][4];

	state_grid[0][0] = state[0];
	state_grid[0][1] = state[1];
	state_grid[0][2] = state[2];
	state_grid[0][3] = state[3];
	state_grid[1][0] = state[4];
	state_grid[1][1] = state[5];
	state_grid[1][2] = state[6];
	state_grid[1][3] = state[7];
	state_grid[2][0] = state[8];
	state_grid[2][1] = state[9];
	state_grid[2][2] = state[10];
	state_grid[2][3] = state[11];
	state_grid[3][0] = state[12];
	state_grid[3][1] = state[13];
	state_grid[3][2] = state[14];
	state_grid[3][3] = state[15];

	addroundkey(state_grid,0,sub,ekey);
	loop_main : for(iteration = 1; iteration < nr; iteration++) {
	#pragma HLS loop_tripcount min=13 max=13
	#pragma HLS pipeline
	   subbytes(sub,shift);
	   shift_row_enc(shift,mix);
	   mixcolumn(mix,round);
	   addroundkey(round,iteration,sub,ekey);
	  }
	  subbytes(sub,shift);
	  shift_row_enc(shift,round);
	  addroundkey(round,nr,result,ekey);

	  cipher[0] = result[0][0];
	  cipher[1] = result[0][1];
	  cipher[2] = result[0][2];
	  cipher[3] = result[0][3];
	  cipher[4] = result[1][0];
	  cipher[5] = result[1][1];
	  cipher[6] = result[1][2];
	  cipher[7] = result[1][3];
	  cipher[8] = result[2][0];
	  cipher[9] = result[2][1];
	  cipher[10] = result[2][2];
	  cipher[11] = result[2][3];
	  cipher[12] = result[3][0];
	  cipher[13] = result[3][1];
	  cipher[14] = result[3][2];
	  cipher[15] = result[3][3];

}

void shift_row_enc(uint8_t state[4][4], uint8_t result[4][4])
{


	result[0][0] = state[0][0];
	result[0][1] = state[0][1];
	result[0][2] = state[0][2];
	result[0][3] = state[0][3];

	result[1][0] = state[1][1];
	result[1][1] = state[1][2];
	result[1][2] = state[1][3];
	result[1][3] = state[1][0];

	result[2][0] = state[2][2];
	result[2][2] = state[2][0];
	result[2][1] = state[2][3];
	result[2][3] = state[2][1];

	result[3][0] = state[3][3];
	result[3][3] = state[3][2];
	result[3][2] = state[3][1];
	result[3][1] = state[3][0];
}

void subbytes(uint8_t state[4][4], uint8_t result[4][4])
{

	uint8_t x, y; //addresses the matrix
	loop_sb1 : for(x=0;x<4;x++){
	#pragma HLS loop_tripcount min=4 max=4
	#pragma HLS unroll factor=4 skip_exit_check
		loop_sb2 : for(y=0;y<4;y++){
		#pragma HLS loop_tripcount min=4 max=4
		#pragma HLS unroll factor=4 skip_exit_check
			result[x][y] = sbox[state[x][y]];
		}//end y loop
	}//end x loop
}

void addroundkey(uint8_t state[4][4], uint8_t iteration,uint8_t result[4][4],uint8_t ekey[240])
{
	  #pragma HLS array_partition variable=ekey complete
	  uint8_t x,y;
	  loop_rk1 :for(x=0;x<4;x++) {
	  #pragma HLS loop_tripcount min=4 max=4
	  #pragma HLS unroll factor=4 skip_exit_check
		  	loop_rk2 : for(y=0;y<4;y++){
			#pragma HLS loop_tripcount min=4 max=4
			#pragma HLS unroll factor=4 skip_exit_check
			  result [y][x] = state [y][x] ^ ekey[iteration * nb * 4 + x * nb + y];
	    }
	  }
}

void mixcolumn(uint8_t state[4][4],uint8_t result[4][4])
{

#pragma HLS inline

	//B1’ = (B1 * 2) XOR (B2 * 3) XOR (B3 * 1) XOR (B4 * 1)
	//B2’ = (B1 * 1) XOR (B2 * 2) XOR (B3 * 3) XOR (B4 * 1)
	//B3’ = (B1 * 1) XOR (B2 * 1) XOR (B3 * 2) XOR (B4 * 3)
	//B4’ = (B1 * 3) XOR (B2 * 1) XOR (B3 * 1) XOR (B4 * 2)
	uint8_t x; // control of the column
	loop_mx1 :for(x=0;x<4;x++){
	#pragma HLS loop_tripcount min=4 max=4
	#pragma HLS unroll factor=4 skip_exit_check

		result[0][x] = (gf2[state[0][x]])^(gf3[state[1][x]])^(state[2][x])^(state[3][x]);
		result[1][x] = (state[0][x])^(gf2[state[1][x]])^(gf3[state[2][x]])^(state[3][x]);
		result[2][x] = (state[0][x])^(state[1][x])^(gf2[state[2][x]])^(gf3[state[3][x]]);
		result[3][x] = (gf3[state[0][x]])^(state[1][x])^(state[2][x])^(gf2[state[3][x]]);
	}
}
