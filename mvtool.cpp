//
//  mvtool.cpp
//
//  Copyright 2023 Marco Oliva. All rights reserved.
//

#include <CLI/CLI.hpp>

// from pfp
#include <vcf.hpp>

// from leviosam
#include <leviosam.hpp>

void
serialize_lifting_curr_sample(vcfbwt::Sample& sample, std::size_t sample_genotype, std::size_t length, std::ofstream& out_stream)
{
    lift::Lift_builder lvs_builder(length);
    
    for(size_t v = 0; v < sample.variations.size(); v++)
    {
        const vcfbwt::Variation& variation = sample.get_variation(v);
        
        const size_t var_genotype = sample.genotypes[v][sample_genotype];
        // Get only variation in the current genotype
        if( var_genotype == 0) { continue; }
        
        int rlen = variation.alt[0].size(); // Length of the reference allele
        int alen = variation.alt[var_genotype].size(); // Length of the alternate allele
        lvs_builder.set(variation.pos, variation.types[var_genotype], rlen, alen );
    }
    
    // Build lifting data_structure
    lift::Lift lift(lvs_builder);
    size_t offset = 0;
    out_stream.write((char *)&offset, sizeof(offset));
    // Serialize the data structure
    lift.serialize(out_stream);
}

int main(int argc, char **argv)
{
    CLI::App app("MVTool");
    
    std::size_t w = 10;
    std::vector<std::string> vcfs_file_names;
    std::vector<std::string> refs_file_names;
    std::size_t max_samples = 0;
    std::string samples_file_name;
    std::string haplotype_string = "1";
    std::string out_prefix;
    std::string tmp_dir;
    
    app.add_option("-v,--vcf", vcfs_file_names, "List of comma ',' separated vcf files. Assuming in genome order!")->allow_extra_args(true)->configurable()->delimiter(',')->required();
    app.add_option("-r,--ref", refs_file_names, "List of comma ',' separated reference files. Assuming in genome order!")->allow_extra_args(true)->configurable()->delimiter(',')->required();
    app.add_option("-w, --window-size", w, "Sliding window size.")->check(CLI::Range(3, 200))->configurable()->required();
    app.add_option("-m, --max", max_samples, "Max number of samples to analyze")->configurable();
    app.add_option("-o,--out-prefix", out_prefix, "Output prefix")->configurable()->required();
    app.add_option("-S, --samples", samples_file_name, "File containing the list of samples to parse")->configurable();
    app.add_option("-H,--haplotype", haplotype_string, "Haplotype: [1,2,12].")->configurable();
    app.add_option("--tmp-dir", tmp_dir, "Temporary files directory.")->check(CLI::ExistingDirectory)->configurable();
    app.allow_windows_style_options();

    CLI11_PARSE(app, argc, argv);
    
    // clean file name vectors
    vcfs_file_names.erase(std::remove_if(vcfs_file_names.begin(), vcfs_file_names.end(),
                                         [] (std::string& s) {return s.size() == 0; } ), vcfs_file_names.end());
    refs_file_names.erase(std::remove_if(refs_file_names.begin(), refs_file_names.end(),
                                         [] (std::string& s) {return s.size() == 0; } ), refs_file_names.end());
    
    // print out configurations
    spdlog::info("Current Configuration:\n{}", app.config_to_str(true,true));
    
    // set tmp file dir
    if (not tmp_dir.empty()) { vcfbwt::TempFile::setDirectory(tmp_dir); }
    
    // parse the VCF
    vcfbwt::VCF vcf(refs_file_names, vcfs_file_names, samples_file_name, max_samples);
    
    // generate indexes
    std::vector<std::string> names;
    std::vector<std::size_t> lengths;
    
    std::string tmp_lifting_name = vcfbwt::TempFile::getName("lifting");
    std::ofstream tmp_lifting_out_stream(tmp_lifting_name);
    
    // reference length
    std::string reference_id = "reference";
    std::size_t reference_length =  vcf.get_reference().size() + w;
    
    names.push_back(reference_id);
    lengths.push_back(reference_length);
    
    for (std::size_t i = 0; i < vcf.size(); i++)
    {
        if (haplotype_string == "1" or haplotype_string == "2")
        {
            std::size_t sample_genotype;
            if (haplotype_string == "1") { sample_genotype = 0; }
            else { sample_genotype = 1; }
            
            vcfbwt::Sample::iterator it(vcf[i], sample_genotype);
            
            // lengths
            std::size_t length = it.length() + w;
            if (i == vcf.size() - 1) { length += w - 1; } // last sample has w dollar primes and a dollar
            names.push_back(vcf[i].id());
            lengths.push_back(length);
            
            // lifting
            serialize_lifting_curr_sample(vcf[i], sample_genotype, length, tmp_lifting_out_stream);
        }
        else
        {
            // first haplotype
            vcfbwt::Sample::iterator it_h1(vcf[i], 0);
    
            // second haplotype
            vcfbwt::Sample::iterator it_h2(vcf[i], 1);
            
            // lengths
            std::size_t length_h1 = it_h1.length() + w;
            std::size_t length_h2 = it_h2.length() + w;
            if (i == vcf.size() - 1) { length_h2 += w - 1; } // last sample has w dollar primes and w dollars
    
            names.push_back("H1_" + vcf[i].id());
            lengths.push_back(length_h1);
    
            names.push_back("H2_" + vcf[i].id());
            lengths.push_back(length_h2);
    
            // lifting
            serialize_lifting_curr_sample(vcf[i], 0, length_h1, tmp_lifting_out_stream);
            serialize_lifting_curr_sample(vcf[i], 1, length_h2, tmp_lifting_out_stream);
        }
    }
    
    // write out lengths
    std::ofstream lengths_out_stream(out_prefix + ".lengths");
    for (std::size_t i = 0; i < names.size(); i++)
    {
        lengths_out_stream << names[i] << " " << lengths[i] << std::endl;
    }
    
    vcfbwt::DiskWrites::update(lengths_out_stream.tellp());
    lengths_out_stream.close();
    
    
    // Writing out lifitng structure
    vcfbwt::DiskWrites::update(tmp_lifting_out_stream.tellp());
    tmp_lifting_out_stream.close();
    
    std::vector<std::size_t> onset(1,0);
    std::size_t u = 0;
    
    for (std::size_t i = 0; i < names.size(); i++)
    {
        u += lengths[i];
        onset.push_back(u);
    }
    ++u;
    
    std::ofstream lifting_out_stream(out_prefix + ".lifting");
    
    // seqidx structure
    sdsl::sd_vector_builder builder(u, onset.size());
    for (auto idx : onset) { builder.set(idx); }
    
    sdsl::sd_vector<> starts(builder);
    sdsl::sd_vector<>::rank_1_type rank1(&starts);
    sdsl::sd_vector<>::select_1_type select1(&starts);
    
    lifting_out_stream.write((char *)&u, sizeof(u));
    lifting_out_stream.write((char *)&w, sizeof(w));
    
    starts.serialize(lifting_out_stream);
    sdsl::serialize(names.size(), lifting_out_stream);
    for(size_t i = 0; i < names.size(); ++i)
    {
        sdsl::serialize(names[i].size(), lifting_out_stream);
        lifting_out_stream.write((char *)names[i].data(), names[i].size());
    }
    
    // total number of contigs
    std::size_t n_contigs = 1;
    lifting_out_stream.write((char *)&n_contigs, sizeof(n_contigs));
    
    // empty lifting for the references
    const size_t len = vcf.get_reference().size() + w;
    
    sdsl::bit_vector ibv(len);
    sdsl::bit_vector dbv(len);
    sdsl::bit_vector sbv(len);
    
    lift::Lift lift(ibv, dbv, sbv);
    
    sdsl::serialize(len, lifting_out_stream);
    lift.serialize(lifting_out_stream);
    
    // merge tmp liftings from vcf
    std::ifstream tmp_lifting_in(tmp_lifting_name);
    lifting_out_stream << tmp_lifting_in.rdbuf();
    tmp_lifting_in.close();
    
    vcfbwt::DiskWrites::update(lifting_out_stream.tellp());
    lifting_out_stream.close();
}