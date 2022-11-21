#include "OpenGrowth.h"

// Let's say that we want to choose a number between 0 and 19, each one having a given probability to be picked (stored in the array probaFile).
// This function returns one of the number between 0 and 19 with the appropriate probability.
int Random(Parameters & parameters, double const probaFile[])
{
    int randomValue=0;
    int choiceHasBeenMade=0;

    // We do this kind of loop because with approximate floating numbers, the sum of all the probabilities may not be exactly 1.0000.
    // If the sum of probabilities is 0.99980000 and the randomNumber is 0.99985000, we will not pick a randomValue and do the loop
    // again. The other reason is that it allows the user to change a probability file and put some numbers to 0.00000000. The relative
    // probabilities between two choices will not change, and the program will still work.
    while ( choiceHasBeenMade==0 ) {
        // We take a random number between 0 and 99999999 and divide by 100000000 to have it between 0.00000000 and 0.99999999.
        long double randomNumber = (parameters.random() % 100000000) / 100000000.0 ;
        long double sum=0;
        for (int j=0 ; j < parameters.fragmentListSize ; j++) {
            if ( (sum<=randomNumber) && (randomNumber<(sum+probaFile[j])) ) { randomValue=j; j=parameters.fragmentListSize-1; choiceHasBeenMade++; }
            sum += probaFile[j];
            }
        }

    return randomValue;
}

