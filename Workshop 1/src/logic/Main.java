package logic;

import java.io.IOException;
import java.util.Scanner;

public class Main {
	public static void main(String[] args) throws IOException{
		int option = 3;
		int n = 0;
		int m = 0;
		String nombre="";
		boolean validInput = false;
		Scanner scanner = new Scanner(System.in);
		do {
			System.out.print("\n--Menu-- \n");
			System.out.print("1. Create database \n"
			                 + "2. Exit\n");
			do {
				try {
					System.out.print("Input a number: ");
					option = Integer.parseInt(scanner.nextLine());
					if (option > 2) {
						System.out.println("That is not a valid number. Try again.");
					} else {
						validInput = true;
					}
				} catch (NumberFormatException e) {
                System.out.println("That is not a valid number. Try again.");
                }
			}while (validInput != true);
			switch(option) {
			case 1:
				 
				System.out.print("Enter the number n of sequences of dataset:(1000 <= n <= 2000000): ");
				do {
					try {
						n = Integer.parseInt(scanner.nextLine());
						if (n < 1000 || n > 2000000 ) {
							System.out.println("That is not a valid number. Try again.");
							validInput = false;
						} else {
							validInput = true;
							break;
						}
					} catch (NumberFormatException e) {
	                System.out.println("That is not a valid number. Try again.");
	                validInput = false;
	                }
				}while (!validInput);
				
				System.out.print("Enter the number m of size of the motiv:(4 <= m <= 10): ");
				do {
					try {
						
						m = Integer.parseInt(scanner.nextLine());
						if (m < 4 || m > 10) {
							System.out.println("That is not a valid number. Try again.");
							validInput = false;
						} else {
							validInput = true;
						}
					} catch (NumberFormatException e) {
	                System.out.println("That is not a valid number. Try again.");
	                validInput = false;
	                }
				}while (!validInput);
				Bioinformatics bio = new Bioinformatics(n, 5, 100, 0.25, 0.25, 0.25, 0.25);
				bio.createDataset();
				nombre = "dataset "+String.valueOf(n)+" "+String.valueOf(m);
				bio.saveDataset(nombre,m);
				String mostFrequentMotif = bio.getMotif(m);
                System.out.println("Most frequent motif of length " + m + ": " + mostFrequentMotif);
				System.out.println("Saved dataset data without taking shannon entropy into account in the .txt called: "+nombre+" with a quantity of:"+bio.getQuantity()+"\n");
				bio.loadDataset(nombre);
				bio.applyEntropy();
				nombre = "new "+nombre;
				bio.saveModifiedDataset(nombre,m);
				mostFrequentMotif = bio.getMotif(m);
                System.out.println("Most frequent motif of length " + m + ": " + mostFrequentMotif);

				System.out.println("Saved dataset data  taking shannon entropy into account in the .txt called:"+nombre+" with a quantity of:"+bio.getQuantity()+"\n");
				break;
			default:
				break;
			}
		}while(option != 2);
    }

}
