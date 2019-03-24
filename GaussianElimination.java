//Scott Ha
//CS 3010
//Assignment 1
package numerical;

import java.io.File;
import java.util.Random;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Stack;
import java.util.Scanner;

public class GaussianElimination {
    // creates a long start, end, and total time variable to time an alg.
    int num;
    float[][] matrix1;
    float[][] Matrix2;
    Scanner scan;
  
    // constructor
    public GaussianElimination() throws FileNotFoundException{
        scan = new Scanner(System.in);
        num = getNum();
        scan.nextLine();
        matrix1 = getMatrix();
        Matrix2 = Arrays.copyOf(matrix1, matrix1.length);
    } 
  
    //Takes user input for number of equations 
    private int getNum() {
        System.out.print("Enter the number of equations: ");
        return scan.nextInt();
    } 
  
    //Takes user choice of input for coefficients for FILE or Command Line 
    private float[][] getMatrix() throws FileNotFoundException {
        matrix1 = new float[num][num+1];
        int choice;
        System.out.print("Enter '1' to input coefficients by command line, or '2' to enter by filename: ");
        choice = scan.nextInt();
        scan.nextLine();
        if(choice == 1) commandLine();
        else fileInput();      
        return matrix1;
    } 
  
    //takes command line input for Coefficients
    private void commandLine(){
        String[] values;
        //Random rand = new Random();
        for (int i = 0; i < num; i++) {
            System.out.print("Enter the " + (num+1) +
                    " coefficients for row " + (i+1) +
                    " (seperate coefficients by space): ");
            values = scan.nextLine().split("\\s+");
            for(int j = 0; j < values.length; j++){
                matrix1[i][j] = Float.parseFloat(values[j]);
            }
        }
    } 
  
    //Takes file input for coefficient matrix
    private void fileInput() throws FileNotFoundException {
        System.out.print("Enter the filename: ");
        File file = new File(scan.nextLine());
        try (Scanner scan = new Scanner(file)) {
            for (int i = 0; i < num && scan.hasNextFloat(); i++) {
                for(int j = 0; j < num+1 && scan.hasNextFloat(); j++){
                    matrix1[i][j] = scan.nextFloat();
                }
            }
            scan.close();
        }
        catch(FileNotFoundException e){}
    }
  
    //Scaled Partial Pivoting for Gaussian Elimination
    private void PartialPivoting(float[][] matrix1){
    	
        float[] weights = getWeights(matrix1);
        int pivotRow;
        float[] results = new float[num];
        Arrays.fill(results, 0);      
        Stack<Integer> pivots = new Stack<>();
        for(int i = 0; i < num; i++){
            // get pivotal equation
            pivotRow = getPivotRow(i,weights,pivots,matrix1);
            gaussElimination(i,pivotRow,weights,pivots,matrix1);
        }
        backSub(matrix1,pivots,results);
        printResults(results);
    }
  
    //Gets weights by interating through matrix
    private float[] getWeights(float[][] matrix1) {
        float[] weights = new float[num];
        float maximum;
        for(int i = 0; i < num; i++) {
            maximum = -1;
            for(int j = 0; j < num+1; j++) {
                maximum = Math.max(maximum, Math.abs(matrix1[i][j]));
            }
            weights[i] = maximum;
        }
        return weights;
    }
  
    //Gets max weight and pivot ratio 
    private int getPivotRow(int i, float[] weights, Stack<Integer> pivots, float[][] matrix){
        int pivotRow;
        float[][] maxRatioPair = { {-1, -1} };
        float ratio;
        for(int j = 0; j < num; j++){
            if(!pivots.contains(j)) {
                if(i+1 == num) {
                    pivots.push(j);
                    break;
                }
                else {
                    ratio = Math.abs(matrix[j][i] / weights[j]);
                    if(ratio > maxRatioPair[0][1]) {
                        maxRatioPair[0][0] = j;
                        maxRatioPair[0][1] = ratio;
                    }
                }
            }
        }
        if(i+1 != num){
            pivotRow = (int) maxRatioPair[0][0];
            pivots.push(pivotRow);  
        }
        else pivotRow = pivots.peek();
        return pivotRow;
    }
  
    //Gaussian Elimination method
    private void gaussElimination(int i, int pivotRow,float[] weights, Stack<Integer> pivots, float[][] matrix1){
        float targetNumber;
        float pivotWeight;
    	float pivotNumber;
        float pivotRate;
        if(i+1 != num) {
            for(int j = 0; j < num; j++) {
                if((!pivots.contains(j)) || matrix1[j][i] == 0) {
                    pivotRate = matrix1[pivotRow][i] / matrix1[j][i];
                    pivotWeight = weights[pivotRow];
                    for(int k = 0; k < num+1; k++) {
                        pivotNumber = matrix1[pivotRow][k];
                        targetNumber = matrix1[j][k];
                        matrix1[j][k] = (pivotNumber - pivotRate * targetNumber) / pivotWeight;
                    }
                    printMatrix(matrix1);
                }
            }
        }
    } 
  
    //Back substitution for Gaussian Method
    private void backSub(float[][] matrix1, Stack<Integer> pivots, float[] results){      
        	int pivotRow;
        	float tempWeight;
        for(int i = num-1; !pivots.empty();i--) {
            pivotRow = pivots.pop();
            tempWeight = matrix1[pivotRow][i];
            results[i] = matrix1[pivotRow][num] / tempWeight;
            for(int j = num-1; j >= 0 ; j--) {
                if(i != j){
                    results[i] -= (matrix1[pivotRow][j] / tempWeight) * results[j];
                }
            }
        }
    } 
    
    //Jacobi Method 
    private void Jacobi(float[][] matrix1){
        // asks user for desired error
        System.out.print("Enter the desired error in decimal form: ");
        float error = scan.nextFloat();
        float totalError = 0;
        float[] results = new float[num];
        float[] previousResults;
        float[] tempResults = new float[num];
        Arrays.fill(results, 0);
        // print intemediate x coefficient matrix
        System.out.println("x^" + 0 + ":");
        printResults(results);
        // 50 iterations if desired error not achieved
        for(int i = 0; i < 50; i++) {
            previousResults = results.clone();
            for(int j = 0; j < num; j++){
                tempResults[j] = matrix1[j][num] / matrix1[j][j];
                // column
                for(int k = num-1; k >= 0; k--){
                    if(k != j) {
                        tempResults[j] -= results[k] * (matrix1[j][k] / matrix1[j][j]);
                    }
                }
            }
            results = tempResults.clone();
            tempResults = new float[num];
            if(i > 0){
                for (int k = 0; k < results.length; k++){
                    totalError += Math.abs(results[k] - previousResults[k]) / 2;
                }
                if(totalError < error){
                    System.out.println("x^" + (i+1) + ":");
                    printResults(results);
                    break;
                }
                totalError = 0;
            }
            System.out.println("x^" + (i+1) + ":");
            printResults(results);
        }
    } 
  

    	//Gauss Siedel Method
    private void Siedel(float [][] matrix1){
        // asks user for desired error
        System.out.print("Enter the desired error in decimal form: ");
        float error = scan.nextFloat();
        float totalError = 0;
        float[] results = new float[num];
        float[] previousResults;
        Arrays.fill(results, 0);
        // print intemediate x coefficient matrix
        System.out.println("x^" + 0 + ":");
        printResults(results);
        // 50 iterations if desired error not achieved
        for(int i = 0; i < 50; i++) {
            // store old results
            previousResults = Arrays.copyOf(results, results.length);
            for(int j = 0; j < num; j++){
                results[j] = matrix1[j][num] / matrix1[j][j];
                for(int k = num-1; k >= 0; k--){
                    if(k != j) {
                        results[j] -= results[k] * (matrix1[j][k] / matrix1[j][j]);
                    }
                }
            }
            // will only calculate error when there's a previous iteration
            if(i > 0){
                for (int k = 0; k < results.length; k++){
                    totalError += Math.abs(results[k] - previousResults[k]) / 2;
                }
                if(totalError < error){
                    System.out.println("x^" + (i+1) + ":");
                    printResults(results);
                    break;
                }
                totalError = 0;
            }
            System.out.println("x^" + (i+1) + ":");
            printResults(results);
        }
    } 
  

    
    //Print Function for results
    public void printResults(float[] results){
        DecimalFormat df = new DecimalFormat("#.###");
        for(int i = 0; i < results.length; i++){
            System.out.println("x" + (i+1) + " = " + df.format(results[i]));
        }
        System.out.println("");
    } 
  
    //Print function for matrix
    public static void printMatrix(float[][] matrix1){
        DecimalFormat df = new DecimalFormat("#.###");
        for (float[] row : matrix1) {
            System.out.print("[ ");
            for (int column = 0; column < row.length; column++) {
                System.out.print(df.format(row[column]) + " ");
            }
            System.out.println("]");
        }
        System.out.println("");
    }     
          
    public static void Matrix(float[][] target, float[][] source){
        for(int i = 0; i < source.length; i++)
            target[i] = source[i].clone();
    }
  
   
    //Main function; includes the 3 functions Partial Pivoting, Gauss Jacobi, Gauss Siedel
    public static void main(String[] args) throws FileNotFoundException {
        GaussianElimination G = new GaussianElimination();
        float[][] Matrix1 = G.matrix1;
        float[][] Matrix2 = new float[Matrix1.length][Matrix1[0].length];
        Matrix(Matrix2, Matrix1);
        GaussianElimination.printMatrix(Matrix2);
      
        System.out.println("Scaled Partial Pivoting Gaussian Elimination:");
        G.PartialPivoting(Matrix2);
        Matrix(Matrix2, Matrix1);
      
        System.out.println("Gauss Jacobi:");
        G.Jacobi(Matrix2);
        Matrix(Matrix2, Matrix1);
      
        System.out.println("Gauss Siedel:");
        G.Siedel(Matrix2);
    } 
  
} 