package logic;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Vector;

public class Bioinformatics {
	
	private double[] weights;
	private int datasetSize;
	private Vector<String> dataset;
	private int minLenght;
	private int maxLenght;
	private long duration;
	private double entropy;
	private int quantity;
	
	public double getQuantity() {
		return quantity;
	}

	public void setQuantity(int quantity) {
		this.quantity = quantity;
	}

	public double getEntropy() {
		return entropy;
	}

	public void setEntropy(double entropy) {
		this.entropy = entropy;
	}

	public Bioinformatics(int datasetSize, int minLenght, int maxLenght,
		                 double weightA, double weightC, double weightG, double weightT) {
		this.datasetSize = datasetSize;
		this.minLenght = minLenght;
		this.maxLenght = maxLenght;
		this.weights = new double[4];
		this.dataset = new Vector<>();
		this.weights[0]= weightA;
		this.weights[1]= weightC;
		this.weights[2]= weightG;
		this.weights[3]= weightT;
		this.duration = 0;
		this.createDataset();
		
	}
	
	public String createSequence(int size) {
		StringBuilder sequence = new StringBuilder();
		Random ran = new Random();
		for (int i = 0; i < size; i++) {
			double value = ran.nextDouble();
			if (value < weights[0]) {
				sequence.append("A");
			} else if (value < weights[0] + weights[1]) { 
				sequence.append("C");
			} else if (value < weights[0] + weights[1] + weights[2]) { 
				sequence.append("G");
			} else {
				sequence.append("T");
			}
		}
		return sequence.toString();
	}
	
	public void createDataset() {
		Random ran = new Random();
		for (int i = 0; i < this.datasetSize; i++) {
			int sequenceSize = this.minLenght + ran.nextInt(this.maxLenght - this.minLenght + 1);
			dataset.add(createSequence(sequenceSize));
		}
	}
	
	public String getMotif(int sizeMotif) {
		long  startTime = System.currentTimeMillis();
		HashMap<String, Integer> motifMap = new HashMap<>();
		String mostMotif = null;
		int maxCount = 0;
		for (String sequence : dataset) {
			for (int i = 0; i <= sequence.length() - sizeMotif; i++) {
				String motif = sequence.substring(i, i + sizeMotif);
				int count = motifMap.getOrDefault(motif, 0) + 1;
				motifMap.put(motif, count);
				if (count > maxCount) {
					maxCount = count;
					mostMotif = motif;
				}
			}
		}
		this.quantity = maxCount;
		long  endTime = System.currentTimeMillis();
		long sequenceDuration = endTime - startTime;
		this.duration += sequenceDuration;
		return mostMotif;
	}
	
	public void saveDataset(String filename, int motifLength) {
		try {
			Path directoryPath = Paths.get("data");
			if (!Files.exists(directoryPath)) {
				Files.createDirectory(directoryPath); 
			}
			File file = new File(directoryPath.toFile(), filename);

			try (FileWriter writer = new FileWriter(file)) {
				for (String sequence : dataset) {
					writer.write(sequence + "\n");
				}
				String mostFrequentMotif = getMotif(motifLength);
				writer.write("\nMost frequent motif of length " + motifLength + ": " + mostFrequentMotif + " with a quantity of:"+this.quantity+"\n");
				writer.write("Time to find motif: " + this.duration + " ms\n");

			} catch (IOException e) {
				System.err.println("Error writing to file: " + e.getMessage());
			}
		} catch (IOException e) {
			System.err.println("Error creating directory: " + e.getMessage());
		}
	}
	
	public void loadDataset(String filename) {
		dataset.clear(); 
		try {
			Path filePath = Paths.get("data", filename);
			Files.lines(filePath).forEach(line -> {
				if (!line.trim().isEmpty() && !line.startsWith("Most frequent motif")) {
					dataset.add(line); 
				}
			});
		} catch (IOException e) {
			System.err.println("Error reading from file: " + e.getMessage());
		}
	}
	
	public char generateRandomBase(char currentBase) {
	    Random random = new Random();
	    char newBase = 'A'; 
	    do {
	        int randomValue = random.nextInt(4);
	        switch (randomValue) {
	            case 0:
	                newBase = 'A';
	                break;
	            case 1:
	                newBase = 'C';
	                break;
	            case 2:
	                newBase = 'G';
	                break;
	            case 3:
	                newBase = 'T';
	                break;
	        }
	    } while (newBase == currentBase); 

	    return newBase;
	}
	
	public String modifySequenceByEntropy(String sequence, double entropy) {
		StringBuilder newSequence = new StringBuilder(sequence);
		for (int i = 0; i < sequence.length(); i++) {
		    if (entropy < 1.5) { 
		        char currentBase = sequence.charAt(i);
		        char newBase = generateRandomBase(currentBase); 
		        newSequence.setCharAt(i, newBase);
		    }
		}
		return newSequence.toString();
	}
	
	public void saveModifiedDataset(String filename , int motifLength) {
		try {
			Path directoryPath = Paths.get("data");
			if (!Files.exists(directoryPath)) {
				Files.createDirectory(directoryPath);
			}
			File file = new File(directoryPath.toFile(), filename);
			try (FileWriter writer = new FileWriter(file)) {
				for (String sequence : dataset) {
					writer.write(sequence + "\n");
				}
				writer.write("Dataset modified by entropy\n");
				writer.write("Entropy: "+this.getEntropy()+"\n");
				String mostFrequentMotif = getMotif(motifLength);
				writer.write("\nMost frequent motif of length " + motifLength + ": " + mostFrequentMotif + " with a quantity of:"+this.quantity+"\n");
				writer.write("Time to find motif: " + this.duration + " ms\n");

			} catch (IOException e) {
				System.err.println("Error writing to file: " + e.getMessage());
			}
		} catch (IOException e) {
			System.err.println("Error creating directory: " + e.getMessage());
		}
	}

	public void applyEntropy() {
		Vector<String> newDataset = new Vector<>();
		for (String sequence : dataset) {
			double entropy = calculateEntropy(sequence);
			String newSequence = modifySequenceByEntropy(sequence, entropy); 
			newDataset.add(newSequence);
		}
		dataset = newDataset; 
	}
	
	public double calculateEntropy(String sequence) {
        Map<Character, Integer> freqMap = new HashMap<>();
        for (char base : sequence.toCharArray()) {
            freqMap.put(base, freqMap.getOrDefault(base, 0) + 1);
        }

        double entropy = 0.0;
        int length = sequence.length();
        for (int count : freqMap.values()) {
            double probability = (double) count / length;
            entropy -= probability * Math.log(probability) / Math.log(2);
        }
        this.setEntropy(entropy);
        return entropy;
    }
}
