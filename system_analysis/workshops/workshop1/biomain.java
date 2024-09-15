package system_analysis.workshops.workshop1;

// author: Angel Diaz
public class biomain {
    public static void main(String[] args) {
        int min_dataset_size = 1000;
        int max_dataset_size = 2000000;
        int min_length = 5;
        int max_length = 100;
        double weight_a = 0.25;
        double weight_c = 0.25;
        double weight_g = 0.25;
        double weight_t = 0.25;
        int motif_size = 4;

        BioInformatics bio = new BioInformatics(min_dataset_size, max_dataset_size, min_length, max_length, weight_a,
                weight_c, weight_g, weight_t, motif_size);

        int action = -1;
        do {
            java.util.Scanner scanner = new java.util.Scanner(System.in);
            System.out.println("Select an action:");
            System.out.println("1: Generate dataset");
            System.out.println("2: Load dataset");
            System.out.println("3: Get motif from dataset");
            System.out.println("4: Filter by Shannon entropy");
            System.out.println("5: load filtered dataset");
            System.out.println("6: Get motif from filtered dataset");
            System.out.println("0: Exit");

            while (!scanner.hasNextInt()) {
                System.out.println("Please enter a valid integer.");
                scanner.next(); // consume the invalid input
            }
            action = scanner.nextInt();

            switch (action) {
                case 0:
                    break;
                case 1:
                    bio.generate_dataset();
                    System.out.println("dataset generated.");
                    break;
                case 2:
                    bio.load_dataset(1);
                    System.out.println("Dataset loaded.");
                    break;
                case 3:
                    bio.get_motif(1);
                    System.out.println("Motif obtained.");
                    break;
                case 4:
                    bio.shannon_entropy_filter();
                    System.out.println("Filtered by Shannon entropy obtained");
                    break;
                case 5:
                    bio.load_dataset(2);
                    System.out.println("Filtered dataset loaded.");
                    break;
                case 6:
                    bio.get_motif(2);
                    System.out.println("filtered Motif obtained.");
                    break;
                default:
                    System.out.println("Invalid action.");
                    break;
            }

        } while (action != 0);

    }

}
