package my;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

import org.biojava.bio.BioException;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequence.IOTools;
import org.biojavax.bio.seq.RichSequenceIterator;
import org.jcvi.jillion.core.qual.QualitySequence;
import org.jcvi.jillion.trace.fastq.FastqFileReader;
import org.jcvi.jillion.trace.fastq.FastqQualityCodec;

import utils.CompareSequence;
import utils.Subject;

public class my_SIQ {
    private static Subject subject;

    public static void main(String[] args) {
        File inputFile = new File(args[0]);
        File refFile = new File(args[1]);
        String leftFlank = args[2];
        String rightFlank = args[3];
        File outputFile = new File(args[4]);
        try(BufferedReader is = new BufferedReader(new FileReader(refFile))) {
            RichSequenceIterator si = IOTools.readFastaDNA(is, null);
            try {
                RichSequence rs = si.nextRichSequence();
                subject = new Subject(rs, leftFlank, rightFlank);
                subject.swapPrimersIfNeeded();
            } catch(BioException e) {
                e.printStackTrace();
            }
        } catch(IOException e) {
            e.printStackTrace();
        }
			
		try (PrintWriter writer = new PrintWriter(new FileOutputStream(outputFile), false)) {
            FastqFileReader.forEach(
                inputFile,
                FastqQualityCodec.SANGER, 
                (id, fastqRecord) -> {
                    QualitySequence quals = fastqRecord.getQualitySequence();
                    String seq = fastqRecord.getNucleotideSequence().toString();
                    
                    CompareSequence cs = new CompareSequence(subject, seq, quals, inputFile.getParentFile().getName(), false, id);
                    cs.setDelinsFilter(true);
                    
                    cs.setAndDetermineCorrectRange(0.08);
                    cs.maskSequenceToHighQualityRemoveSingleRange();
                    cs.setAllowJump(false);
                    
                    boolean leftCorrect = false;
                    boolean rightCorrect = false;
                    //at this point has to be true because of earlier check
                    if(cs.isMasked()) {
                        cs.setCurrentAlias(null, inputFile.getName());
                        cs.determineFlankPositions(true);
                        leftCorrect = cs.isCorrectPositionLeft();
                        rightCorrect = cs.isCorrectPositionRight();
                        //only correctly found ones
                        if(leftCorrect && rightCorrect){
                            int rpos1 = cs.leftFlank.getSubjectEnd();
                            int rpos2 = cs.subjectObject.subject.indexOf(cs.rightFlank, rpos1);
                            int qpos1 = cs.leftFlank.getQueryEnd();
                            int qpos2 = cs.query.indexOf(cs.rightFlank, qpos1);
                            
                            writer.println(id+"\t"+rpos1+"\t"+rpos2+"\t"+qpos1+"\t"+qpos2);
                        }
                    }
                }
            );
        }
        catch(IOException e) {
            e.printStackTrace();
        }
    }
}
