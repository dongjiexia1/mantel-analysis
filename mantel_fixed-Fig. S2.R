library(readxl)
library(linkET)
library(tidyverse)

# ============================================
# 1. 读取两个独立的文件 (图1)
# ============================================
spec <- read_excel("C:/Users/17652/Desktop/相关性分析（夏）/图1/root_architecture.xlsx") %>%
  column_to_rownames("Sample") %>%
  mutate(across(everything(), as.numeric))

env  <- read_excel("C:/Users/17652/Desktop/相关性分析（夏）/图1/1retu.xls") %>%
  column_to_rownames("Sample") %>%
  mutate(across(everything(), as.numeric))

# ============================================
# 2. 从分组信息.xlsx 自动读取分组（超出spec范围的列自动剔除）
# ============================================
group_info <- read_excel("C:/Users/17652/Desktop/相关性分析（夏）/图1/分组信息.xlsx")
spec_group <- split(group_info$Indicator, group_info$Group)
spec_group <- lapply(spec_group, function(x) intersect(x, colnames(spec)))
spec_group <- spec_group[lengths(spec_group) > 0]
cat("Loaded groups:", paste(names(spec_group), collapse = ", "), "\n")

# ============================================
# 3. Mantel test (p 值矫正)
# ============================================
mantel <- mantel_test(spec, env,
                       spec_select = spec_group) %>%
    rename(p_raw = p) %>%
    mutate(p_adj = p.adjust(p_raw, method = "BH"),
           rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                    labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
           pd = cut(p_adj, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")),
           r_sign = ifelse(r >= 0, "+", "−"))
corr <- correlate(env, method = "spearman")
corr_p_raw <- corr$p   # 保存原始 p
corr$p <- corr$p %>%
  as.matrix() %>%
  p.adjust(method = "BH") %>%
  matrix(nrow = nrow(corr$p), dimnames = dimnames(corr$p)) %>%
  as.data.frame()

# ============================================
# 4. 画图
# ============================================
p <- qcorrplot(corr, type = "upper") +
  geom_square() +
  geom_text(
    data = {
      r_mat <- corr$r
      p_mat <- corr$p
      vnames <- colnames(r_mat)
      expand.grid(row = vnames, col = vnames) %>%
        filter(as.character(row) != as.character(col)) %>%
        rowwise() %>%
        mutate(
          p_value = p_mat[as.character(row), as.character(col)],
          sig = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01  ~ "**",
            p_value < 0.05  ~ "*",
            TRUE            ~ ""
          )
        ) %>%
        ungroup() %>%
        filter(match(row, vnames) < match(col, vnames))
    },
    aes(x = col, y = row, label = sig),
    size = 4, color = "black", fontface = "bold",
    inherit.aes = FALSE
  ) +
  geom_couple(aes(colour = pd, size = rd, linetype = r_sign), data = mantel,
              curvature = 0.15) +
  scale_linetype_manual(values = c("+" = "solid", "−" = "dashed"),
                        name = "Mantel's r sign") +
  scale_fill_gradient2(low  = "#2166AC", mid  = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman's r") +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size   = guide_legend(title = "Mantel's r",
                               override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p",
                               override.aes = list(size = 3), order = 1),
         fill   = guide_colorbar(title = "Spearman's r", order = 3))

# 输出
ggsave("C:/Users/17652/Desktop/相关性分析（夏）/图1/Mantel_plot_Fig1.png", plot = p,
       width = 14, height = 12, dpi = 300, units = "in")
ggsave("C:/Users/17652/Desktop/相关性分析（夏）/图1/Mantel_plot_Fig1.svg", plot = p,
       width = 14, height = 12, dpi = 300, units = "in")

# ============================================
# 5. 导出表格
# ============================================
corr_r <- as.data.frame(corr$r) %>%
  rownames_to_column("env1") %>%
  pivot_longer(-env1, names_to = "env2", values_to = "Spearman_r")

corr_p <- as.data.frame(corr$p) %>%
  rownames_to_column("env1") %>%
  pivot_longer(-env1, names_to = "env2", values_to = "p_adj")

corr_p_raw_long <- as.data.frame(corr_p_raw) %>%
  rownames_to_column("env1") %>%
  pivot_longer(-env1, names_to = "env2", values_to = "p_raw")

spearman_table <- left_join(corr_r, corr_p_raw_long, by = c("env1", "env2")) %>%
  left_join(corr_p, by = c("env1", "env2")) %>%
  filter(match(env1, colnames(corr$r)) < match(env2, colnames(corr$r))) %>%
  mutate(Significance = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01  ~ "**",
    p_adj < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))

mantel_table <- mantel %>% select(spec, env, r, p_raw, p_adj, rd, pd)

library(openxlsx)
write.xlsx(list(Spearman = spearman_table, Mantel = mantel_table),
           "C:/Users/17652/Desktop/相关性分析（夏）/图1/Mantel_results_Fig1.xlsx")

cat("图1 Done!\n")
