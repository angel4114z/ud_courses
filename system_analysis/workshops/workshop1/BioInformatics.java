package system_analysis.workshops.workshop1;

import java.util.Random;
import java.util.Vector;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;


public class BioInformatics {

    private double[] weights = {0, 0, 0, 0}; // probability of each nucleotide
    private int dataset_size; // number of sequences in the dataset
    private Vector<String> dataset = new Vector<String>(); // dataset of sequences
    private int min_length; // minimum length of a sequence
    private int max_length; // maximum length of a sequence
    private int motif_size; // size of the motif
    private HashMap<String, Integer> motif = new HashMap<String, Integer>();
    private ExecutorService executor;

    public BioInformatics(int min_dataset_size,
                          int max_dataset_size,
                          int min_length,
                          int max_length,
                          double weight_a,
                          double weight_c,
                          double weight_g,
                          double weight_t,
                          int motif_size
                        ) {
        this.dataset_size = new Random().nextInt(max_dataset_size - min_dataset_size) + min_dataset_size;
        this.weights[0] = weight_a;
        this.weights[1] = weight_c;
        this.weights[2] = weight_g;
        this.weights[3] = weight_t;
        this.min_length = min_length;
        this.max_length = max_length;
        this.motif_size = motif_size;
    }

    public void generate_dataset() {

        executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

        for (int i = 0; i < this.dataset_size; i++) {
            executor.submit(() -> generate_secuence());
        }

        executor.shutdown();
        shutdownExecutor();

        save_dataset(1);
        
    }

    private void generate_secuence(){

        Random rand = new Random();
        int length = rand.nextInt(this.max_length - this.min_length) + this.min_length;
            String sequence = "";
            double weight[] = {this.weights[0], this.weights[1] + this.weights[0] , this.weights[2] + this.weights[1] + this.weights[0]  , this.weights[3]};

            for (int j = 0; j < length; j++) {
                double r = rand.nextDouble();
                if (r < weight[0]) {
                    sequence += "A";
                } else if (r < weight[1]) {
                    sequence += "C";
                } else if (r < weight[2]) {
                    sequence += "G";
                } else {
                    sequence += "T";
                }
            }
            
            synchronized(this.dataset){
                this.dataset.add(sequence);
            }
    
        }

    private void save_dataset(int q) {
        try {

            String filename = "";
            if(q == 1){
                filename = "dataset.txt";
            }else if(q == 2){
                filename = "filtered_dataset.txt";
            }

            FileWriter writer = new FileWriter(filename);
            
            
            for (String sequence : this.dataset) {
                writer.write(sequence + "\n");
            }
            writer.close();
            System.out.println("Dataset saved successfully.");
        } catch (IOException e) {
            System.out.println("Error saving dataset: " + e.getMessage());
        }
    }

    public void load_dataset(int q) {
        try {
            this.dataset.clear();

            String filename = "";
            if(q == 1){
                filename = "dataset.txt";
            }else if(q == 2){
                filename = "filtered_dataset.txt";
            }

            java.io.BufferedReader reader = new java.io.BufferedReader(new java.io.FileReader(filename));
            String line;
            while ((line = reader.readLine()) != null) {
                this.dataset.add(line);
            }
            reader.close();
            System.out.println("Dataset loaded successfully.");
        } catch (IOException e) {
            System.out.println("Error loading dataset: " + e.getMessage());
        }
    }


    public void get_motif(int q) {

        dataset_size = this.dataset.size();

        System.out.println("dataset size: " + dataset_size);

        long startTime = System.currentTimeMillis();
        
        executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

        for(int i = 0 ; i < this.dataset.size(); i++){
            final int index = i;
            executor.submit(() -> find_motif(index));
        }

        executor.shutdown();
        shutdownExecutor();
        

        long endTime = System.currentTimeMillis();
        long totalTime = endTime - startTime;
        long seconds = totalTime / 1000;
        System.out.println("Total time taken: " + totalTime + " milliseconds\n" + seconds + " seconds");

        save_motif(totalTime, q);
        
    }

    private void save_motif(long totalTime, int q) {

        int max_motif = 0;
        String max_motif_key = "";

        long seconds = totalTime / 1000;

        for(String key : motif.keySet()){
            if(motif.get(key) > max_motif){
                max_motif = motif.get(key);
                max_motif_key = key;
            }
        }
        System.out.println("The most repeated motif is: " + max_motif_key + " with " + max_motif + " repetitions");

        try {
            String filename = "";
            if(q == 1){
                filename = "motif_counts.txt";
            }else if(q == 2){
                filename = "filtered_motif_counts.txt";
            }
            FileWriter writer = new FileWriter(filename);
            writer.write("dataset size: " + this.dataset_size + "\n");
            writer.write("motif size: " + this.motif_size + "\n");
            writer.write("most repeated motif: " + max_motif_key + " with " + max_motif + " repetitions\n");
            writer.write("probability of each nucleotide: A=" + this.weights[0] + ", C=" + this.weights[1] + ", G=" + this.weights[2] + ", T=" + this.weights[3] + "\n");
            writer.write("total time taken: " + totalTime + " milliseconds\n " + seconds + " seconds\n");

            for (String key : motif.keySet()) {
                int count = motif.get(key);
                writer.write(key + ": " + count + "\n");
            }
            writer.close();
            System.out.println("Motif counts saved successfully.");
        } catch (IOException e) {
            System.out.println("Error saving motif counts: " + e.getMessage());
        }

    }

    private void find_motif(int i){

        motif.clear();
        String temp_secuence = this.dataset.get(i);
        int temp_size = temp_secuence.length();

        System.out.println("secuence size: " + temp_size);

        for(int j = 0 ; j <= temp_size - this.motif_size; j += this.motif_size){


            String motif_candidate = temp_secuence.substring(j, j + this.motif_size);
            if(motif.containsKey(motif_candidate)){
                motif.put(motif_candidate, motif.get(motif_candidate)+1);
            }else{
                motif.put(motif_candidate, 1);
            }
        }


    }

    public void shannon_entropy_filter() {
        dataset_size = this.dataset.size();

        executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        long startTime = System.currentTimeMillis();
        
        for(int i = 0 ; i < this.dataset.size(); i++){
            final int index = i;
            executor.submit(() -> shannon_entropy(index));
        }

        executor.shutdown();
        shutdownExecutor();

        long endTime = System.currentTimeMillis();
        long totalTime = endTime - startTime;
        long seconds = totalTime / 1000;
        System.out.println("Total time taken: " + totalTime + " milliseconds\n" + seconds + " seconds");

        save_dataset(2);

        
    }

    private void shannon_entropy(int i){
        String temp_secuence = this.dataset.get(i);
        String replace_secuence = "";
        int temp_size = temp_secuence.length();

        for(int j = 0 ; j < temp_size - this.motif_size; j += this.motif_size){

            if (j + this.motif_size < temp_size) {

            String temp_motif = temp_secuence.substring(j, j + this.motif_size);

            
            int[] count = new int[4];
            
            for (char nucleotide : temp_motif.toCharArray()) {
                switch (nucleotide) {
                    case 'A': count[0]++; break;
                    case 'C': count[1]++; break;
                    case 'G': count[2]++; break;
                    case 'T': count[3]++; break;
                }
            }
            double[] probabilities = new double[4];
            
            probabilities[0] = (double) count[0] / this.motif_size;
            probabilities[1] = (double) count[1] / this.motif_size;
            probabilities[2] = (double) count[2] / this.motif_size;
            probabilities[3] = (double) count[3] / this.motif_size;
            
            
            double entropy = 0;
            for (double p : probabilities) {
                if (p > 0.0) {
                    entropy -= (p * (Math.log(p) / Math.log(2)));
                }
            }
            
            
            if(entropy > 1.6){
                replace_secuence += temp_motif;
            }
            
        }

        this.dataset.set(i, replace_secuence);
        }
    }
        

    private void shutdownExecutor() {
        executor.shutdown();
        try {
            if (!executor.awaitTermination(300, TimeUnit.SECONDS)) {
                executor.shutdownNow();
            }
        } catch (InterruptedException e) {
            executor.shutdownNow();
        }
    }
}