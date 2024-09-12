package system_analysis.workshops.workshop1;

import java.util.Random;
import java.util.Vector;
import java.util.concurrent.CyclicBarrier;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;


public class BioInformatics {

    private double[] weights = {0, 0, 0, 0}; // probability of each nucleotide
    private int dataset_size; // number of sequences in the dataset
    private Vector<String> dataset = new Vector<String>(); // dataset of sequences
    private int min_length; // minimum length of a sequence
    private int max_length; // maximum length of a sequence
    private int motif_size; // size of the motif
    private CyclicBarrier barrera;

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
        this.weights[1] = weight_c + weight_a;
        this.weights[2] = weight_g + this.weights[1];
        this.weights[3] = weight_t + this.weights[2];
        this.min_length = min_length;
        this.max_length = max_length;
        this.motif_size = motif_size;
    }

    public void generate_dataset() {

        System.out.println("Generating dataset with " + this.dataset_size + " sequences");

        barrera = new CyclicBarrier(this.dataset_size);

        for (int i = 0; i < this.dataset_size; i++) {
            Thread t = new Thread(() -> generate_secuence());
            t.start();
        }
        
        try {
            System.out.println("Waiting for threads to finish");
            barrera.await();
            System.out.println("All threads finished");
        } catch (InterruptedException e) {
            System.out.println("interruped \n Error: " + e.getMessage());
        } catch (BrokenBarrierException e) {
            System.out.println("barrier \n Error: " + e.getMessage());
        }
        
    }

    public void generate_secuence(){

        System.out.println(" generate Thread " + Thread.currentThread().getName() + " started");

        Random rand = new Random();
        int length = rand.nextInt(this.max_length - this.min_length) + this.min_length;
            String sequence = "";
            for (int j = 0; j < length; j++) {
                double r = rand.nextDouble();
                if (r < this.weights[0]) {
                    sequence += "A";
                } else if (r < this.weights[1]) {
                    sequence += "C";
                } else if (r < this.weights[2]) {
                    sequence += "G";
                } else {
                    sequence += "T";
                }
            }

            System.out.println("generate Thread " + Thread.currentThread().getName() + " finished");

            this.dataset.add(sequence);
            try {
                barrera.await();
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (BrokenBarrierException e) {
                e.printStackTrace();
            }
            
    }

    public void save_dataset() {
        try {
            FileWriter writer = new FileWriter("dataset.txt");
            for (String sequence : this.dataset) {
                writer.write(sequence + "\n");
            }
            writer.close();
            System.out.println("Dataset saved successfully.");
        } catch (IOException e) {
            System.out.println("Error saving dataset: " + e.getMessage());
        }
    }


    public void get_motif() {
        HashMap<String, Integer> motif = new HashMap<String, Integer>();
        
        for(int i = 0 ; i < this.dataset.size(); i++){
            String temp_secuence = this.dataset.get(i);
            int temp_size = temp_secuence.length();
            for(int j = 0 ; j < temp_size - this.motif_size; j++){
                String motif_candidate = temp_secuence.substring(j, j + this.motif_size);
                if(motif.containsKey(motif_candidate)){
                    motif.put(motif_candidate, motif.get(motif_candidate)+1);
                }else{
                    motif.put(motif_candidate, 1);
                }
            }
        }
        int max_motif = 0;
        String max_motif_key = "";

        for(String key : motif.keySet()){
            if(motif.get(key) > max_motif){
                max_motif = motif.get(key);
                max_motif_key = key;
            }
        }
        System.out.println("The most repeated motif is: " + max_motif_key + " with " + max_motif + " repetitions");

        try {
            FileWriter writer = new FileWriter("motif_counts.txt");
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




    public static void main(String[] args) {
        BioInformatics bio = new BioInformatics(5,100, 1000, 2000000, 0.25, 0.25, 0.25, 0.25, 4);
        bio.generate_dataset();

        bio.save_dataset();
        //bio.get_motif();

        System.out.println("Program finished");


        
    }
}