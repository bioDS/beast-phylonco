package phylonco.beast.evolution.alignment;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.evolution.alignment.Sequence;
import beastfx.app.inputeditor.AlignmentImporter;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Imports a read-count NEXUS file into a {@link ReadCountAlignment}, so read-count data can be added
 * as a partition in BEAUti.
 *
 * <p>Reads the {@code datatype=readCount} characters block that phylonco-lphy's {@code ReadCountNexus}
 * writes: one matrix row per cell, sites separated by {@code ,} and the four base counts by {@code :}.
 * That is the same encoding {@link ReadCountAlignment} takes in a sequence value, so the row is handed
 * across unchanged.</p>
 *
 * <pre>
 * begin characters;   dimensions nchar=2;
 *     format datatype=readCount;
 *     matrix
 *         cell1 1:0:0:11,0:17:0:12;
 *         cell2 7:0:0:26,0:12:0:8;
 * end;
 * </pre>
 *
 * <p>Only the counts matrix is read; the taxa block is redundant with it. mpileup is deliberately not
 * handled here — its parsing lives in lphy-base, and phylonco-beast does not depend on lphy.</p>
 */
@Description("Imports read-count NEXUS files (datatype=readCount) as a read-count alignment")
public class ReadCountNexusImporter implements AlignmentImporter {

    private static final String DATATYPE_MARKER = "datatype=readcount";

    @Override
    public String[] getFileExtensions() {
        // no leading dot: AlignmentImporter compares against the text after the last '.'
        return new String[]{"nexus", "nex"};
    }

    /**
     * Claims a file only if it actually declares the read-count datatype, so ordinary NEXUS
     * alignments are left to the standard parser.
     */
    @Override
    public boolean canHandleFile(File file) {
        if (!AlignmentImporter.super.canHandleFile(file)) {
            return false;
        }
        try {
            for (String line : Files.readAllLines(file.toPath(), StandardCharsets.UTF_8)) {
                if (normalise(line).contains(DATATYPE_MARKER)) {
                    return true;
                }
            }
        } catch (IOException e) {
            return false;
        }
        return false;
    }

    @Override
    public List<BEASTInterface> loadFile(File file) {
        List<String> lines;
        try {
            lines = Files.readAllLines(file.toPath(), StandardCharsets.UTF_8);
        } catch (IOException e) {
            throw new IllegalArgumentException("could not read " + file + ": " + e.getMessage(), e);
        }

        ReadCountAlignment alignment = new ReadCountAlignment();
        boolean sawDataType = false;
        boolean inMatrix = false;
        int rows = 0;

        for (String raw : lines) {
            String line = stripComments(raw).trim();
            if (line.isEmpty()) {
                continue;
            }
            String lower = normalise(line);

            if (lower.contains(DATATYPE_MARKER)) {
                sawDataType = true;
            }
            if (!inMatrix) {
                // "matrix" may sit alone on its line or lead the first data row
                int at = lower.indexOf("matrix");
                if (at >= 0) {
                    inMatrix = true;
                    line = line.substring(at + "matrix".length()).trim();
                    if (line.isEmpty()) {
                        continue;
                    }
                } else {
                    continue;
                }
            }
            if (lower.startsWith("end;") || line.equals(";")) {
                break;
            }

            String[] parts = line.split("\\s+", 2);
            if (parts.length != 2) {
                throw new IllegalArgumentException(
                        "malformed read-count matrix row in " + file.getName() + ": '" + raw.trim() + "'");
            }
            String counts = parts[1].trim();
            if (counts.endsWith(";")) {
                counts = counts.substring(0, counts.length() - 1);
            }
            alignment.setInputValue("sequence", new Sequence(parts[0], counts));
            rows++;
        }

        if (!sawDataType) {
            throw new IllegalArgumentException(
                    file.getName() + " is not a read-count NEXUS file (no 'datatype=readCount')");
        }
        if (rows == 0) {
            throw new IllegalArgumentException("no read-count matrix rows found in " + file.getName());
        }

        alignment.setID(idFrom(file));
        alignment.initAndValidate();

        List<BEASTInterface> result = new ArrayList<>();
        result.add(alignment);
        return result;
    }

    private static String normalise(String line) {
        return stripComments(line).trim().toLowerCase(Locale.ROOT).replace(" ", "");
    }

    /** Drops NEXUS [ ... ] comments. */
    private static String stripComments(String line) {
        return line.replaceAll("\\[[^\\]]*\\]", "");
    }

    private static String idFrom(File file) {
        String name = file.getName();
        int dot = name.lastIndexOf('.');
        return dot > 0 ? name.substring(0, dot) : name;
    }
}
